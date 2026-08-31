function kine = compute_kinematics(fish_points, fps, min_freq)
% COMPUTE_KINEMATICS  FFT-based kinematic analysis on transformed fish midlines.
%
%   (See original docstring, unchanged from prior version. This file
%   patches two local helpers, fill_nan, dominant_freq, so that
%   degenerate/all-NaN input produces NaN output instead of a silently
%   fabricated number, fixes fft_interp's duplicated-endpoint bug,
%   converts Y to the world-frame lateral displacement with the
%   rigid-body recoil removed, measures tail_TBF from the raw
%   camera-space tail displacement (see CHANGE NOTE in section 4), and
%   fits the wavelength with an explicit complex-offset traveling-wave
%   model in spatial_wavelength (see its CHANGE NOTE).)

    nFish = numel(fish_points);
    N_OUT = 200;

    if isscalar(fps), fps = repmat(fps, nFish, 1); end

    kine(nFish) = struct();

    for fi = 1:nFish

        X_raw  = fish_points(fi).X;
        Y_raw  = fish_points(fi).Y;
        has_z  = isfield(fish_points(fi), 'Z') && ~isempty(fish_points(fi).Z) ...
                 && isfield(fish_points(fi), 'has_z') && fish_points(fi).has_z;
        Z_raw  = [];
        if has_z, Z_raw = fish_points(fi).Z; end

        nFrames = size(X_raw, 1);
        fs = fps(fi);
        s_norm = linspace(0, 1, N_OUT);

        % ---- Early check: warn loudly if there's no usable data at all ----
        n_frames_with_data = sum(~any(isnan(X_raw), 2) & ~any(isnan(Y_raw), 2));
        if n_frames_with_data == 0
            warning(['compute_kinematics: %s has ZERO frames with complete X/Y data. ' ...
                     'All outputs for this animal will be NaN (not a fabricated ' ...
                     'placeholder), check transform_fish output / tracking coverage.'], ...
                     fish_points(fi).name);
        end

        % ----------------------------------------------------------------
        % 1.  FFT spatial interpolation: nPoints -> N_OUT per frame
        % ----------------------------------------------------------------
        X_interp = NaN(nFrames, N_OUT);
        Y_interp = NaN(nFrames, N_OUT);
        Z_interp = NaN(nFrames, N_OUT);
        Yc_raw   = NaN(size(Y_raw));   % cleaned raw lateral displacement

        % transform_fish detrends every frame by subtracting its fitted
        % midline (intercept a, slope -tan(theta) through the middle
        % points). That line ITSELF oscillates at the beat frequency, so
        % the raw Y carries a known contamination: the fitted line's
        % value rotated into the body frame. From the transform_fish
        % roundtrip identity Y*bl = (ycam - a + tan(theta)*xcam)*cos(theta)
        % (exact), the fitted line's camera-frame value is
        % (a - tan(theta)*xcam), and its body-frame value is that times
        % cos(theta)/bl. It is added back HERE, at the raw points, before
        % interpolation: interpolating the clean signal is exact for
        % periodic waves, whereas interpolating the contamination and
        % correcting afterwards leaves seam-ringing errors. The result is
        % the WORLD-frame lateral displacement, up to a per-frame constant
        % removed by the spatial-mean subtraction below.
        tp = [];
        if isfield(fish_points(fi), 'transform_params') && ~isempty(fish_points(fi).transform_params)
            tp = fish_points(fi).transform_params;
        end
        has_tp = ~isempty(tp) && all(isfield(tp, {'a','theta','bl','x1'}));
        has_sf = has_tp && all(isfield(tp, {'sign_flip'}));

        for f = 1:nFrames
            x = X_raw(f, :);
            y = Y_raw(f, :);
            z = [];
            if has_z, z = Z_raw(f, :); end

            if any(isnan(x)) || any(isnan(y)), continue; end
            if has_z && any(isnan(z)), continue; end

            if has_tp && all(isfinite([tp(f).a tp(f).theta tp(f).bl tp(f).x1]))
                sf_f = 1;
                if has_sf && isfinite(tp(f).sign_flip), sf_f = tp(f).sign_flip; end
                xcam = (sf_f .* x .* tp(f).bl + tp(f).x1) .* cos(tp(f).theta) ...
                       + y .* tp(f).bl .* sin(tp(f).theta);
                y = y + (tp(f).a - tan(tp(f).theta) .* xcam) .* abs(cos(tp(f).theta)) ./ tp(f).bl;
            end
            Yc_raw(f, :) = y;

            x_tail  = x(end);
            x_query = linspace(0, x_tail, N_OUT);

            X_interp(f, :) = x_query;
            Y_interp(f, :) = fft_interp(y, N_OUT);
            if has_z
                Z_interp(f, :) = fft_interp(z, N_OUT);
            end
        end

        % Remove each frame's spatial mean (the s-independent rigid-body
        % recoil, body-axis sway, tank circling). It adds the same complex
        % offset to every station's phase at the beat frequency and flattens
        % the phase gradient; subtracting it leaves the traveling wave
        % intact (a full-period wave has exactly zero spatial mean, so
        % nothing is lost). Done AFTER interpolation so the mean is
        % estimated on the fine grid, not the few raw points.
        Y_interp = Y_interp - mean(Y_interp, 2, 'omitnan');

        % ----------------------------------------------------------------
        % 2.  Lateral (Y) amplitude envelope
        % ----------------------------------------------------------------
        [amp_mean, amp_std, headAmp, tailAmp, minAmp, minAmpLoc, maxAmp, maxAmpLoc] = ...
            amplitude_stats(Y_interp, s_norm);

        headTailAmpRatio = headAmp / tailAmp;

        % ----------------------------------------------------------------
        % 3.  Dorso-ventral (Z) amplitude (3-D only)
        % ----------------------------------------------------------------
        ampZ_mean = []; ampZ_std = []; headAmpZ = NaN; tailAmpZ = NaN;
        minAmpZ = NaN; minAmpZLoc = NaN; maxAmpZ = NaN; maxAmpZLoc = NaN;
        if has_z
            [ampZ_mean, ampZ_std, headAmpZ, tailAmpZ, minAmpZ, minAmpZLoc, maxAmpZ, maxAmpZLoc] = ...
                amplitude_stats(Z_interp, s_norm);
        end

        % ----------------------------------------------------------------
        % 4.  Beat frequencies: temporal FFT on first (head) and last
        %     (tail) point of the CLEANED lateral displacement.
        %     CHANGE NOTE (bug fix): the raw transform output Y is the
        %     distance from the per-frame fitted midline, and that line
        %     absorbs most of the tail's wave (the line's value at the
        %     tail x rides with the tail beat), so dominant_freq on
        %     Y_raw(:,end) used to return a half-frequency line-tracking
        %     artifact instead of the tail's true lateral beat (a real
        %     shark's tail beat measured 0.92 Hz from Y_raw vs 1.83 Hz
        %     from the cleaned signal AND from the raw camera
        %     coordinates). Yc_raw (section 1) is the world-frame lateral
        %     displacement, so its temporal spectra are the true beats.
        % ----------------------------------------------------------------
        head_Y = fill_nan(Yc_raw(:, 1));
        % CHANGE NOTE (head_TBF trend removal): Yc_raw is the WORLD-frame
        % lateral displacement, so the head's trace also carries the
        % rigid-body translation along the swimming path, a slow trend.
        % dominant_freq's mean removal alone turns a straight-path
        % translation into a sawtooth whose 1/f spectrum buries the beat
        % (a straight-path synthetic at a 30-degree heading measured the
        % min_freq floor instead of its true 2 Hz beat). Remove the
        % best-fit linear trend first: the translation is exactly linear
        % for a straight path, while genuine oscillations (the beat,
        % whole-body sway) survive unchanged.
        if ~all(isnan(head_Y))
            head_Y = detrend(head_Y, 1);
        end
        tail_Y = fill_nan(Yc_raw(:, end));
        % CHANGE NOTE (tail_TBF source): neither the cleaned signal nor
        % the raw camera-space tail y alone is reliable across gait types.
        % On a 5-point benthic walking shark the per-frame line fit
        % (a, theta) oscillates at the STEP frequency (~0.6 Hz), and
        % adding its extrapolated value back at the tail injects that
        % wobble into the cleaned tail signal, burying the real tail beat
        % under a low-frequency artifact (a validated walking trial
        % measured 0.59 Hz from the cleaned signal vs 2.80 Hz from manual
        % tail-beat peak counts). The raw camera-space tail y, in turn,
        % is dominated by whole-body lateral sway on swimmers (a
        % validated swim trial measured 0.61 Hz from it vs a true
        % 1.83 Hz). The tail's camera y RELATIVE TO THE PER-FRAME
        % CENTROID of all tracked points removes the common-mode
        % translation/sway that contaminates the raw signal while
        % carrying none of the midline fit's wobble, it measured the
        % true beat in both validated trials (swim 1.83 Hz, walk
        % 2.95 Hz). The head keeps the cleaned signal: its camera y is
        % dominated by whole-body translation, which the midline-relative
        % Y_raw removes by construction. Pre-transformed data has no
        % camera frame and falls back to Yc_raw, which is then the
        % normalized lateral displacement itself.
        if isfield(fish_points(fi), 'points') ...
           && ~isempty(fish_points(fi).points) ...
           && size(fish_points(fi).points,2) >= 1 ...
           && size(fish_points(fi).points,3) >= 2
            cam_y  = squeeze(fish_points(fi).points(:, :, 2));   % [nF x nPts]
            tail_Y = fill_nan(cam_y(:, end) - mean(cam_y, 2, 'omitnan'));
        end
        [head_TBF, head_freq, head_power] = dominant_freq(head_Y, fs, min_freq);
        [tail_TBF, tail_freq, tail_power] = dominant_freq(tail_Y, fs, min_freq);

        headZ_TBF = NaN;  tailZ_TBF = NaN;
        if has_z
            % CHANGE NOTE (bug fix): Z is normalized so the HEAD station is
            % identically zero (Z = (z - z_head)/bl), so Z_raw(:,1) can never
            % carry a temporal signal and dominant_freq on it always returned
            % NaN. Measure the head-region vertical beat at the station just
            % behind the head instead.
            [headZ_TBF] = dominant_freq(fill_nan(Z_raw(:, min(2, size(Z_raw,2)))), fs, min_freq);
            [tailZ_TBF] = dominant_freq(fill_nan(Z_raw(:, end)), fs, min_freq);
        end

        % ----------------------------------------------------------------
        % 4b. Spline (interpolated-midline) frequency: comparison value
        %     alongside head/tail TBF. Dominant temporal frequency at EACH
        %     interpolated midline station; the summary is the median over
        %     stations that actually oscillate (>= 10% of the max station
        %     amplitude), so near-head stations with tiny/noisy signals
        %     don't skew it. Wave speed itself still uses tail_TBF (below).
        % ----------------------------------------------------------------
        station_freqs = NaN(1, N_OUT);
        station_amps  = NaN(1, N_OUT);
        for k = 1:N_OUT
            y = Y_interp(:, k);
            if all(isnan(y)), continue; end
            station_freqs(k) = dominant_freq(fill_nan(y), fs, min_freq);
            station_amps(k)  = std(y, 'omitnan');
        end
        spline_freq = NaN;
        ok = isfinite(station_freqs) & isfinite(station_amps) & (station_amps > 0);
        if any(ok)
            keep = station_amps >= 0.1 * max(station_amps(ok));
            f_ok = station_freqs(ok & keep);
            if isempty(f_ok), f_ok = station_freqs(ok); end
            spline_freq = median(f_ok);
        end

        % ----------------------------------------------------------------
        % 5.  Propulsive wavelength (two-stage complex-amplitude +
        %     phase-gradient traveling-wave fit, see spatial_wavelength
        %     below for why the old envelope-FFT method was replaced: it
        %     always returned 1.005 BL for every animal, why the
        %     intermediate phase-gradient version was biased for
        %     non-integer wavelengths, and why the final estimate is a
        %     phase gradient on offset-corrected amplitudes)
        % ----------------------------------------------------------------
        f_dom = tail_TBF;
        if isnan(f_dom), f_dom = head_TBF; end
        slin = slinear_contamination(fish_points(fi), nFrames, fs, f_dom);
        [wavelength, wave_sf, wave_pow] = spatial_wavelength(Y_interp, s_norm, fs, f_dom, slin);
        % Wave speed intentionally stays TBF-based: head_TBF / tail_TBF /
        % spline_freq are all exported alongside it for comparison.
        wave_speed = wavelength * tail_TBF;   % BL/s, wave speed = wavelength x TBF

        % ----------------------------------------------------------------
        % 6.  Curvature, XY plane
        % ----------------------------------------------------------------
        [curv_mean, curv_std, maxCurv, maxCurvLoc] = ...
            curvature_stats(X_interp, Y_interp, [], s_norm, nFrames, N_OUT);

        % ----------------------------------------------------------------
        % 7.  Curvature, 3-D (if Z available)
        % ----------------------------------------------------------------
        curv3d_mean = []; curv3d_std = []; maxCurv3D = NaN; maxCurv3DLoc = NaN;
        if has_z
            [curv3d_mean, curv3d_std, maxCurv3D, maxCurv3DLoc] = ...
                curvature_stats(X_interp, Y_interp, Z_interp, s_norm, nFrames, N_OUT);
        end

        % ----------------------------------------------------------------
        % 8.  Store
        % ----------------------------------------------------------------
        kine(fi).name             = fish_points(fi).name;
        kine(fi).X_interp         = X_interp;
        kine(fi).Y_interp         = Y_interp;
        kine(fi).Z_interp         = Z_interp;
        kine(fi).s_norm           = s_norm;
        kine(fi).n_frames_with_data = n_frames_with_data;   % NEW: sanity-check field

        kine(fi).amp_mean         = amp_mean;
        kine(fi).amp_std          = amp_std;
        kine(fi).headAmp          = headAmp;
        kine(fi).tailAmp          = tailAmp;
        kine(fi).headTailAmpRatio = headTailAmpRatio;
        kine(fi).minAmp           = minAmp;
        kine(fi).minAmpLoc        = minAmpLoc;
        kine(fi).maxAmp           = maxAmp;
        kine(fi).maxAmpLoc        = maxAmpLoc;

        kine(fi).ampZ_mean        = ampZ_mean;
        kine(fi).ampZ_std         = ampZ_std;
        kine(fi).headAmpZ         = headAmpZ;
        kine(fi).tailAmpZ         = tailAmpZ;
        kine(fi).minAmpZ          = minAmpZ;
        kine(fi).minAmpZLoc       = minAmpZLoc;
        kine(fi).maxAmpZ          = maxAmpZ;
        kine(fi).maxAmpZLoc       = maxAmpZLoc;

        kine(fi).head_TBF         = head_TBF;
        kine(fi).tail_TBF         = tail_TBF;
        kine(fi).spline_freq_Hz   = spline_freq;   % interpolated-midline dominant freq
        kine(fi).head_fft_freq    = head_freq;
        kine(fi).head_fft_power   = head_power;
        kine(fi).tail_fft_freq    = tail_freq;
        kine(fi).tail_fft_power   = tail_power;
        kine(fi).headZ_TBF        = headZ_TBF;
        kine(fi).tailZ_TBF        = tailZ_TBF;

        kine(fi).wavelength       = wavelength;
        kine(fi).wave_spatial_freq = wave_sf;
        kine(fi).wave_power       = wave_pow;
        kine(fi).wave_speed_BL_s  = wave_speed;   % NEW: wavelength x tail_TBF

        kine(fi).curv_mean        = curv_mean;
        kine(fi).curv_std         = curv_std;
        kine(fi).maxCurv          = maxCurv;
        kine(fi).maxCurvLoc       = maxCurvLoc;

        kine(fi).curv3d_mean      = curv3d_mean;
        kine(fi).curv3d_std       = curv3d_std;
        kine(fi).maxCurv3D        = maxCurv3D;
        kine(fi).maxCurv3DLoc     = maxCurv3DLoc;

        fprintf('%s | %d/%d frames usable | head TBF=%s tail TBF=%s spline=%s wavelength=%s maxCurv=%s\n', ...
                fish_points(fi).name, n_frames_with_data, nFrames, ...
                fmt_val(head_TBF), fmt_val(tail_TBF), fmt_val(spline_freq), ...
                fmt_val(wavelength), fmt_val(maxCurv));
    end
end


% =========================================================================
%  LOCAL HELPERS
% =========================================================================

function s = fmt_val(v)
% Print NaN plainly as "NaN" rather than a misleading "NaN Hz"-style number.
    if isnan(v), s = 'NaN'; else, s = sprintf('%.4f', v); end
end


function y_out = fft_interp(y_in, N_out)
% Zero-pad FFT interpolation from numel(y_in) to N_out points.
%
% CHANGE NOTE (bug fix): y_in is sampled on s = linspace(0,1,N_in), so its
% first and last samples sit at the same physical location s = 0/1 and the
% sequence's true period is N_in-1 samples. The old code FFT'd all N_in
% samples, making the implied period N_in/(N_in-1) too long, a periodic
% signal like sin(2*pi*s) then lands 5% off-bin, leaks across the whole
% spectrum, and interpolates with up to ~30% error near the seam. Dropping
% the duplicated endpoint puts the period back at exactly N_in-1 samples
% (exact interpolation for periodic band-limited signals), and the raw
% endpoint sample is appended for the s=1 station.
    N_in  = numel(y_in);
    if N_in >= N_out
        y_out = interp1(linspace(0,1,N_in), y_in, linspace(0,1,N_out), 'spline');
        return
    end
    yp    = y_in(1:end-1);          % unique samples (drop duplicate at s=1)
    N     = numel(yp);
    M     = N_out - 1;              % interpolate to unique output stations
    Y     = fft(yp);
    half  = floor(N / 2);
    Y_pad = [Y(1:half+1), zeros(1, M - N), Y(half+2:end)];
    yq    = real(ifft(Y_pad)) * (M / N);
    y_out = [yq, y_in(end)];        % s=1 station: keep the raw sample
end


function [amp_mean, amp_std, headAmp, tailAmp, minAmp, minLoc, maxAmp, maxLoc] = ...
         amplitude_stats(D_interp, s_norm)
% Compute amplitude statistics from an interpolated dimension matrix.
    half_amp = abs(D_interp);
    amp_mean = mean(half_amp, 1, 'omitnan');
    amp_std  = std(half_amp, 0, 1, 'omitnan');

    % CHANGE: use 'omitnan' explicitly (was implicit before, a single NaN
    % in the head/tail window used to make headAmp/tailAmp NaN even when
    % most of that window had real data; now it only goes NaN if the WHOLE
    % window is empty of data).
    headAmp = mean(amp_mean(s_norm <= 0.05), 'omitnan');
    tailAmp = mean(amp_mean(s_norm >= 0.95), 'omitnan');

    if all(isnan(amp_mean))
        minAmp = NaN; minLoc = NaN; maxAmp = NaN; maxLoc = NaN;
    else
        [minAmp, mi] = min(amp_mean);  minLoc = s_norm(mi);
        [maxAmp, ma] = max(amp_mean);  maxLoc = s_norm(ma);
    end
end


function [curv_mean, curv_std, maxCurv, maxCurvLoc] = ...
         curvature_stats(X_i, Y_i, Z_i, s_norm, nFrames, N_OUT)
% Per-frame 3-point curvature, averaged across frames.  Z_i may be [].
    use3d    = ~isempty(Z_i);
    curv_all = NaN(nFrames, N_OUT);
    lag      = max(1, round(N_OUT / 40));

    for f = 1:nFrames
        x = X_i(f,:);  y = Y_i(f,:);
        if use3d, z = Z_i(f,:); else, z = zeros(size(x)); end
        if any(isnan(x)) || any(isnan(y)), continue; end
        if use3d && any(isnan(z)), continue; end

        curv_row = NaN(1, N_OUT);
        for k = lag+1 : N_OUT-lag
            x1=x(k-lag); y1=y(k-lag); z1=z(k-lag);
            x2=x(k);     y2=y(k);     z2=z(k);
            x3=x(k+lag); y3=y(k+lag); z3=z(k+lag);

            A = sqrt((x2-x1)^2+(y2-y1)^2+(z2-z1)^2);
            B = sqrt((x3-x2)^2+(y3-y2)^2+(z3-z2)^2);
            C = sqrt((x3-x1)^2+(y3-y1)^2+(z3-z1)^2);
            s = (A+B+C)/2;
            denom = 4*sqrt(max(s*(s-A)*(s-B)*(s-C), 0));
            if denom > 0
                curv_row(k) = 1 / ((A*B*C)/denom);
            end
        end
        curv_all(f,:) = curv_row;
    end

    curv_mean = mean(curv_all, 1, 'omitnan');
    curv_std  = std(curv_all, 0, 1, 'omitnan');

    % CHANGE NOTE: MATLAB's max() on an all-NaN vector returns NaN with
    % index 1 rather than erroring; the ORIGINAL bug used this silently
    % to return curv_mean(1) as if it were a real max. Explicitly check
    % first and propagate NaN/NaN instead of a fake location.
    if all(isnan(curv_mean))
        maxCurv = NaN; maxCurvLoc = NaN;
    else
        [maxCurv, ci] = max(curv_mean);
        maxCurvLoc    = s_norm(ci);
    end
end


function [f_dom, freqs, power] = dominant_freq(y, fs, min_freq)
% CHANGE NOTE (bug fix): previously, if y was entirely NaN, fill_nan()
% (below) silently substituted an all-ZERO vector, and an FFT of a flat
% zero signal always has zero power everywhere; MATLAB's max() on an
% all-zero/all-equal vector returns index 1 by convention instead of
% erroring, so the OLD code always returned the same fake frequency
% (whatever the lowest FFT bin >= min_freq happened to be) with no
% relationship to real data. This version explicitly checks for a
% degenerate signal (all-NaN, or effectively zero variance) FIRST and
% returns NaN instead of a fabricated value.
    if all(isnan(y)) || (max(y) - min(y)) < eps
        f_dom = NaN; freqs = []; power = [];
        return;
    end
    N     = length(y);
    Y     = fft(y - mean(y));
    power = (2/N) * abs(Y(1:floor(N/2)+1)).^2;
    freqs = fs * (0:floor(N/2)) / N;
    valid = freqs >= min_freq;
    if ~any(valid), f_dom = NaN; return; end
    power_valid = power(valid);
    freqs_valid = freqs(valid);
    [~, idx]  = max(power_valid);
    % CHANGE NOTE (leakage gate): when the real dominant frequency lies
    % below min_freq, the valid region holds only that signal's sidelobes
    % (or pure FFT noise for an integer-cycle tone), and the old code
    % returned the strongest of those as a fabricated "frequency", a
    % synthetic 2 Hz beat measured 46 Hz with min_freq = 3. Reject when
    % the valid peak is negligible next to the global peak, or when it is
    % not a local maximum of the FULL spectrum (a stronger neighbor on
    % the invalid side means the valid bin is just a sidelobe of the
    % sub-min_freq peak, not a beat).
    if power_valid(idx) < 1e-8 * max(power)
        f_dom = NaN; return;
    end
    k_peak = idx - 1 + find(freqs >= min_freq, 1);
    if (k_peak > 1 && power(k_peak-1) > power(k_peak)) || ...
       (k_peak < numel(power) && power(k_peak+1) > power(k_peak))
        f_dom = NaN; return;
    end
    f_dom     = freqs_valid(idx);
end


function [wavelength, sf, power] = spatial_wavelength(Y_interp, s_norm, fs, f_dom, slin)
% Wavelength of the propulsive wave via a COMPLEX-AMPLITUDE FIT.
%
% CHANGE NOTE (complete rewrite): the previous version took the FFT of the
% MEAN AMPLITUDE ENVELOPE, mean(|Y|, t). An amplitude envelope is monotonic
% (small at the head, large at the tail) and contains NO phase information,
% so its FFT always peaked at the lowest non-zero spatial-frequency bin,
% every animal on every trial got wavelength = 1/(199/200) = 1.005025 BL,
% identical to six decimals. That value was an artifact of the bin spacing,
% not a measurement.
%
% In a traveling wave each body station oscillates at the beat frequency
% with a phase that lags progressively toward the tail: phi(s) = -k*s with
% k = 2*pi/wavelength. Each station's complex oscillation amplitude G(s)
% at the dominant beat frequency (f_dom) is measured from its temporal
% Fourier coefficient, and the model
%
%       G(s) = c + B*exp(-i*k*s) + d*col_s
%
% is fitted: c is the spatially uniform offset (the per-frame spatial-mean
% removal leaves an offset for any wave with a non-integer number of
% wavelengths per body, its spatial mean oscillates at the beat
% frequency), B the traveling component, and d*col_s the s-linear
% contamination of the camera projection (see the 6th fix below) whose
% SHAPE col_s = -slin*s is known from the transform's line-fit wobble,
% only its complex scale d is free. The fit is linear in the complex
% unknowns (c, B [, d]) for fixed k, so k is scanned over the plausible
% propulsive range (wavelength 0.4-4 BL) and refined with fminbnd.
%
% CHANGE NOTE (4th estimator): the previous version measured each
% station's PHASE and fitted a line through unwrapped phases. A uniform
% complex offset added to every station's Fourier coefficient does not
% shift all phases equally, it perturbs each phase by an amount that
% oscillates along the body, so the phase-gradient fit came out
% systematically biased for non-integer wavelengths (a synthetic
% lambda=0.8 measured 0.820, +2.5% high). Fitting the offset explicitly
% removes the bias: the same synthetic now measures 0.7996, and lambda=1
% improves from 1.0015 to 1.0003.
%
% CHANGE NOTE (5th estimator, two-stage): the complex fit above is
% biased by amplitude MODULATION along the body: a modulated envelope
% leaks into the fitted c and B and pulls k (a synthetic with
% A(s) = A*(1 + 0.5 sin(2 pi s)) measured 0.988 for a true 1.0). The
% final wavelength therefore comes from a second stage: a line fit
% through the unwrapped phases (amplitude^2-weighted), with the fitted
% offset c subtracted first ONLY when the envelope ratio (max/min
% station amplitude over the kept stations) is <= 4, for strongly
% modulated real animals cB(1) is envelope-corrupted and the raw phases
% are used instead. Modulation changes only the amplitudes, not the
% phases, so the phase gradient is envelope-robust either way. This
% combined estimator is exact on the constant-amplitude synthetics
% (1.0000 / 0.8000), nearly unbiased on the modulated one (0.992), and
% matches the manually validated value on a real swim trial (1.2 vs the
% offset-removed stage's 2.6).
%
% CHANGE NOTE (6th fix, s-linear contamination): the correction above
% makes each frame's corrected signal exactly ycam*cos(theta)/bl, but
% the TRUE body lateral is y_b = cos(h)*yc - sin(h)*(xc - cx) + const.
% Expanding around the fitted theta = 2*pi - h - delta(t) (the line-fit
% slope wobbles with the wave), the difference is
%   yc*cos(theta)/bl - y_b/L = (cos(theta)*cos(h) - 1)*y_b/L
%                              + cos(theta)*sin(h)*(-s_true) + const(t).
% The first term only rescales the wave amplitude (its f0 part is
% cos^2(h)*y_b, no phase change), but the second term's f0 coefficient
% is CT*sin(h)*(-s), LINEAR in station s, where CT is the f0
% coefficient of cos(theta(t)). For a 30-degree fish CT ~ 0.005 and the
% contaminant reaches ~48% of the wave amplitude at the tail, dragging
% the phase gradient from -360 to -384 deg/BL (lambda 1.0 measured 0.94).
% The contamination's SHAPE is therefore known, col_s = -CT*sin(h)*s
% with CT measured from the per-frame transform angles and h implied by
% their circular mean, and the complex model is extended with this FIXED
% column (complex scale d fitted): G(s) = c + B*exp(-i*k*s) + d*col_s.
% A free d*s column was tried first and overfit the modulated-envelope
% synthetic (0.955 for a true 1.0); the fixed shape does not, because on
% the heading-0 synthetics CT vanishes identically (cos(theta) wobbles
% only at 2*f0 there) and the column self-nulls, the exact 1.0000/0.8000
% measurements are unaffected, and the tilted fish now measures 1.0000.
%
% AFFINE-CONTAMINATION CORRECTION (2nd fix, kept): Y comes from
% transform_fish, which detrends every frame by fitting a midline
% (intercept a, slope b) through the middle points and subtracting it. But
% that fitted line ITSELF oscillates at the beat frequency; it absorbs a
% significant part of the traveling wave, so each station's Fourier
% amplitude at f_dom is offset by a KNOWN contamination, which pulls every
% station's phase toward a common direction and systematically steepens
% the fitted gradient (a synthetic wave with a true gradient of 450
% deg/BL measured 540 deg/BL -> lambda 0.67 instead of 0.80).
%
% 3rd fix: EXACT correction, applied BEFORE interpolation: the
% contamination is the fitted midline's own value in the rotated body
% frame, (a - tan(theta)*x_cam)*cos(theta)/bl. This is exact: the
% transform_fish roundtrip identity Y*bl = (ycam - a + tan(theta)*xcam)*
% cos(theta) holds to machine precision, so adding the line's value back
% reconstructs ycam*cos(theta)/bl up to a per-frame constant. The
% correction is applied to the RAW points in the interpolation step
% (section 1), interpolating the clean signal is exact for periodic
% waves, whereas interpolating the contamination and correcting
% afterwards leaves seam-ringing errors. Y_interp therefore arrives here
% as the world-frame lateral displacement with each frame's spatial mean
% (the s-independent recoil) already removed.
%
% INPUTS
%   Y_interp - [nFrames x N_OUT] lateral displacement in BL at each body
%              station (from the interpolation step; NaNs allowed)
%   s_norm   - [1 x N_OUT] normalized body position of each station (0..1)
%   fs       - sampling rate (Hz)
%   f_dom    - dominant beat frequency (Hz) at which amplitudes are
%              measured (typically the tail TBF)
%   slin     - (optional) coefficient of the known s-linear contamination
%              CT*sin(h) from slinear_contamination; NaN when unavailable
%
% OUTPUTS
%   wavelength - BL per cycle. NaN if there is no beat frequency, too few
%                usable stations, no oscillation energy, an incoherent
%                fit (standing wave / fin-driven swimmer), or a
%                negligible traveling component (whole body oscillating
%                in phase).
%   sf         - spatial frequency, 1/wavelength in cycles per BL
%   power      - R^2 of the complex-amplitude regression (estimate
%                quality, 0..1). Empty when the fit is rejected, so an
%                empty value always accompanies a NaN wavelength.

    wavelength = NaN; sf = []; power = [];

    nF = size(Y_interp, 1);
    N  = size(Y_interp, 2);
    % Data window: rows with any interpolated data. CSVs often pad the
    % trial with long runs of all-empty rows; the per-station good-frame
    % requirement must be relative to the tracked window, not the padded
    % row total (otherwise a fully-tracked trial in a padded file fails
    % the gate for no reason).
    nF_eff = sum(any(~isnan(Y_interp), 2));
    if isnan(f_dom) || f_dom <= 0 || N < 3 || nF_eff == 0
        return;
    end

    t = (0:nF-1)' / fs;
    w = 2*pi*f_dom;

    % Per-station complex oscillation amplitude at f_dom, evaluated
    % directly (exact even when f_dom falls between FFT bins).
    G = NaN(1, N);
    for k = 1:N
        y = Y_interp(:, k);
        good = ~isnan(y);
        n_good = sum(good);
        if n_good < max(8, nF_eff/5), continue; end
        yg = y(good) - mean(y(good));   % remove DC
        tg = t(good);
        G(k) = (2/n_good) * sum(yg .* exp(-1i*w*tg));
    end

    ok = isfinite(G);
    if sum(ok) < 3, return; end

    s_ok = s_norm(ok); G_ok = G(ok); A_ok = abs(G_ok);

    % No real oscillation anywhere (e.g. a standing wave whose spatial
    % mean removal already zeroed the signal), NaN beats a fabricated
    % wavelength.
    if max(A_ok) < 1e-12, return; end

    % Restrict to stations that actually oscillate (>= 10% of the max
    % station amplitude) so near-head stations with a tiny signal don't
    % wreck the fit with noise.
    keep = A_ok >= 0.1 * max(A_ok);
    if sum(keep) < 3
        keep = true(size(A_ok));   % too few strong stations, use all
    end
    s_ok = s_ok(keep); G_ok = G_ok(keep); A_ok = A_ok(keep);

    % Fit G(s) = c + B*exp(-i*k*s) [+ d*col_s] with amplitude^2 weights.
    % For fixed k the model is linear in the complex unknowns (c, B [,d]);
    % the scan finds the global optimum to scan resolution and fminbnd
    % refines it. col_s is the FIXED s-linear contamination column
    % -slin*s (see the 6th fix header note): its shape is known from the
    % transform's line-fit wobble, so only its complex scale d is free;
    % this keeps the fit from absorbing envelope modulation the way a
    % free d*s column did.
    W     = A_ok(:) .^ 2;  W = W / sum(W);
    s_col = s_ok(:); G_col = G_ok(:);
    col_s = [];
    if nargin >= 5 && isfinite(slin) && abs(slin) > 1e-12
        col_s = -slin .* s_col;
    end
    rss_fun = @(k) complex_fit_rss(k, s_col, G_col, W, col_s);
    lam_scan = linspace(0.4, 4, 2001);
    rss_scan = arrayfun(@(l) rss_fun(2*pi/l), lam_scan);
    [~, i0]  = min(rss_scan);
    k_lo = 2*pi / lam_scan(min(numel(lam_scan), i0+2));
    k_hi = 2*pi / lam_scan(max(1, i0-2));
    k_opt = fminbnd(rss_fun, k_lo, k_hi);
    [rss_opt, cB] = complex_fit_rss(k_opt, s_col, G_col, W, col_s);
    B_opt = cB(2);

    tss = sum(W .* abs(G_col - sum(W .* G_col)).^2);
    r2  = 1 - rss_opt / max(tss, eps);

    % No real traveling wave if the fit is incoherent (r2 < 0.6: the
    % complex amplitudes do not follow a traveling wave, e.g. a standing
    % wave or a fin-driven swimmer whose body midline does not carry a
    % traveling wave) or the traveling component is negligible (< 5% of
    % the strongest station: the whole body oscillating in phase is not
    % a traveling wave). NaN beats a nonsense number.
    if r2 < 0.6 || abs(B_opt) < 0.05 * max(A_ok)
        return;
    end

    % Two-stage estimate: measure the PHASE GRADIENT of the complex
    % amplitudes (amplitude^2-weighted line fit through the unwrapped
    % phases). A uniform complex offset added to every station's Fourier
    % coefficient perturbs each phase by an amount that oscillates along
    % the body, so the offset fitted by the complex model above is
    % subtracted first, BUT only when the envelope is flat enough that
    % the complex fit's offset estimate is trustworthy. A strongly
    % modulated envelope (real animals: near-zero head amplitude growing
    % to the tail) leaks into the fitted c and B and corrupts cB(1), so
    % subtracting it distorts the phases (a validated swim trial measured
    % 2.6 BL vs a manual ~1.2); there the raw phases are used instead;
    % the amplitude^2 weights already de-emphasize the low-amplitude
    % head stations where the offset's phase perturbation is largest.
    % In the flat-envelope branch the fitted s-linear contamination
    % d*col_s is also removed (same corruption caveat applies in the raw
    % branch, where the s-linear contamination is left in place).
    % Modulation changes
    % only the amplitudes, not the phases, so the phase gradient itself
    % is envelope-robust (a synthetic with A(s) = A*(1 + 0.5 sin(2 pi s))
    % measured 0.992 for a true 1.0, where the complex fit alone gave
    % 0.988).
    sw        = sqrt(W);
    Ew        = [sw, sw .* s_col];
    env_ratio = max(A_ok) / max(min(A_ok), eps);
    if env_ratio <= 4
        Gc = G_col - cB(1);
        if ~isempty(col_s), Gc = Gc - cB(3).*col_s; end
        ph = atan2d(-imag(Gc), real(Gc));
    else
        ph = atan2d(-imag(G_col), real(G_col));
    end
    ph_uw = unwrap(ph * pi/180) * 180/pi;           % degrees, continuous
    beta  = Ew \ (sw .* ph_uw);                     % weighted LS: [intercept; slope]
    if abs(beta(2)) > 1e-6
        wavelength = 360 / abs(beta(2));            % BL per cycle
    else
        wavelength = 2*pi / abs(k_opt);             % flat phase, keep complex fit
    end
    sf         = 1 / wavelength;      % cycles per BL
    power      = max(min(r2, 1), 0);
end

function [rss, cB] = complex_fit_rss(k, s_col, G_col, W, col_s)
% Weighted complex least-squares fit of G(s) = c + B*exp(-i*k*s)
% [+ d*col_s] at a fixed spatial frequency k, with amplitude^2 weights W.
% col_s (optional) is the FIXED s-linear contamination column; when
% empty the model is the plain two-term fit and cB has 2 entries.
    if nargin >= 5 && ~isempty(col_s)
        E = [ones(numel(s_col),1), exp(-1i*k*s_col), col_s];
    else
        E = [ones(numel(s_col),1), exp(-1i*k*s_col)];
    end
    sw = sqrt(W);
    cB = (E .* sw) \ (G_col .* sw);
    r  = G_col - E*cB;
    rss = real(sum(W .* abs(r).^2));
end


function slin = slinear_contamination(fish_points, nFrames, fs, f_dom)
% Coefficient CT*sin(h) of the s-linear camera-projection contamination
% (see the 6th fix header note): the f0 coefficient of cos(theta(t))
% times sin of the heading implied by the mean fitted angle
% (theta ~ 2*pi - h). NaN when the transform params are unavailable.
    slin = NaN;
    if isnan(f_dom) || f_dom <= 0, return; end
    if ~isfield(fish_points, 'transform_params') ...
            || isempty(fish_points.transform_params), return; end
    tp = fish_points.transform_params;
    if ~all(isfield(tp, {'theta'})), return; end
    th = [tp.theta].';
    if numel(th) ~= nFrames, return; end
    good = isfinite(th);
    if sum(good) < 8, return; end
    t = (0:nFrames-1)' / fs;
    w = 2*pi*f_dom;
    CT = (2/sum(good)) * sum(cos(th(good)) .* exp(-1i*w*t(good)));
    thm = angle(sum(exp(1i*th(good))));       % circular mean of theta
    slin = CT * sin(2*pi - thm);
end


function y = fill_nan(y)
% CHANGE NOTE (bug fix): previously, an all-NaN column was silently
% replaced with an all-ZERO vector ("y(:) = 0"). That flat-zero signal
% would then flow into dominant_freq() and produce a fake, deterministic
% "frequency" with no connection to real data (see dominant_freq above).
% Now an all-NaN column stays NaN, and dominant_freq()'s own guard
% catches it and returns NaN cleanly instead.
    t    = (1:length(y))';
    good = ~isnan(y);
    if sum(good) < 2
        y(:) = NaN;
        return;
    end
    y(~good) = interp1(t(good), y(good), t(~good), 'linear', 'extrap');
end