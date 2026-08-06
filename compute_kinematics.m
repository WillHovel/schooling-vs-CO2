function kine = compute_kinematics(fish_points, fps, min_freq)
% COMPUTE_KINEMATICS  FFT-based kinematic analysis on transformed fish midlines.
%
%   (See original docstring — unchanged from prior version. This file
%   patches three local helpers — fill_nan, dominant_freq, spatial_wavelength
%   — so that degenerate/all-NaN input produces NaN output instead of a
%   silently fabricated number. See CHANGE NOTE below each helper.)

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

        [nFrames, nPoints] = size(X_raw);
        fs = fps(fi);
        s_norm = linspace(0, 1, N_OUT);

        % ---- Early check: warn loudly if there's no usable data at all ----
        n_frames_with_data = sum(~any(isnan(X_raw), 2) & ~any(isnan(Y_raw), 2));
        if n_frames_with_data == 0
            warning(['compute_kinematics: %s has ZERO frames with complete X/Y data. ' ...
                     'All outputs for this animal will be NaN (not a fabricated ' ...
                     'placeholder) — check transform_fish output / tracking coverage.'], ...
                     fish_points(fi).name);
        end

        % ----------------------------------------------------------------
        % 1.  FFT spatial interpolation — nPoints -> N_OUT per frame
        % ----------------------------------------------------------------
        X_interp = NaN(nFrames, N_OUT);
        Y_interp = NaN(nFrames, N_OUT);
        Z_interp = NaN(nFrames, N_OUT);

        for f = 1:nFrames
            x = X_raw(f, :);
            y = Y_raw(f, :);
            z = [];
            if has_z, z = Z_raw(f, :); end

            if any(isnan(x)) || any(isnan(y)), continue; end
            if has_z && any(isnan(z)), continue; end

            x_tail  = x(end);
            x_query = linspace(0, x_tail, N_OUT);

            X_interp(f, :) = x_query;
            Y_interp(f, :) = fft_interp(y, N_OUT);
            if has_z
                Z_interp(f, :) = fft_interp(z, N_OUT);
            end
        end

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
        % 4.  Beat frequencies — temporal FFT on first (head) and last (tail) point
        % ----------------------------------------------------------------
        head_Y = fill_nan(Y_raw(:, 1));
        tail_Y = fill_nan(Y_raw(:, end));
        [head_TBF, head_freq, head_power] = dominant_freq(head_Y, fs, min_freq);
        [tail_TBF, tail_freq, tail_power] = dominant_freq(tail_Y, fs, min_freq);

        headZ_TBF = NaN;  tailZ_TBF = NaN;
        if has_z
            [headZ_TBF] = dominant_freq(fill_nan(Z_raw(:,1)),   fs, min_freq);
            [tailZ_TBF] = dominant_freq(fill_nan(Z_raw(:,end)), fs, min_freq);
        end

        % ----------------------------------------------------------------
        % 5.  Propulsive wavelength
        % ----------------------------------------------------------------
        mean_Y_profile = mean(abs(Y_interp), 1, 'omitnan');
        [wavelength, wave_sf, wave_pow] = spatial_wavelength(mean_Y_profile, s_norm);

        % ----------------------------------------------------------------
        % 6.  Curvature — XY plane
        % ----------------------------------------------------------------
        [curv_mean, curv_std, maxCurv, maxCurvLoc] = ...
            curvature_stats(X_interp, Y_interp, [], s_norm, nFrames, N_OUT);

        % ----------------------------------------------------------------
        % 7.  Curvature — 3-D (if Z available)
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
        kine(fi).n_frames_with_data = n_frames_with_data;   % NEW — sanity-check field

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
        kine(fi).head_fft_freq    = head_freq;
        kine(fi).head_fft_power   = head_power;
        kine(fi).tail_fft_freq    = tail_freq;
        kine(fi).tail_fft_power   = tail_power;
        kine(fi).headZ_TBF        = headZ_TBF;
        kine(fi).tailZ_TBF        = tailZ_TBF;

        kine(fi).wavelength       = wavelength;
        kine(fi).wave_spatial_freq = wave_sf;
        kine(fi).wave_power       = wave_pow;

        kine(fi).curv_mean        = curv_mean;
        kine(fi).curv_std         = curv_std;
        kine(fi).maxCurv          = maxCurv;
        kine(fi).maxCurvLoc       = maxCurvLoc;

        kine(fi).curv3d_mean      = curv3d_mean;
        kine(fi).curv3d_std       = curv3d_std;
        kine(fi).maxCurv3D        = maxCurv3D;
        kine(fi).maxCurv3DLoc     = maxCurv3DLoc;

        fprintf('%s | %d/%d frames usable | head TBF=%s tail TBF=%s wavelength=%s maxCurv=%s\n', ...
                fish_points(fi).name, n_frames_with_data, nFrames, ...
                fmt_val(head_TBF), fmt_val(tail_TBF), fmt_val(wavelength), fmt_val(maxCurv));
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
    N_in  = numel(y_in);
    if N_in >= N_out
        y_out = interp1(linspace(0,1,N_in), y_in, linspace(0,1,N_out), 'spline');
        return
    end
    Y     = fft(y_in);
    half  = floor(N_in / 2);
    Y_pad = [Y(1:half+1), zeros(1, N_out - N_in), Y(half+2:end)];
    y_out = real(ifft(Y_pad)) * (N_out / N_in);
end


function [amp_mean, amp_std, headAmp, tailAmp, minAmp, minLoc, maxAmp, maxLoc] = ...
         amplitude_stats(D_interp, s_norm)
% Compute amplitude statistics from an interpolated dimension matrix.
    half_amp = abs(D_interp);
    amp_mean = mean(half_amp, 1, 'omitnan');
    amp_std  = std(half_amp, 0, 1, 'omitnan');

    % CHANGE: use 'omitnan' explicitly (was implicit before — a single NaN
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
    % index 1 rather than erroring — the ORIGINAL bug used this silently
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
% zero signal always has zero power everywhere — MATLAB's max() on an
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
    f_dom     = freqs_valid(idx);
end


function [wavelength, sf, power] = spatial_wavelength(y_profile, s_norm)
% CHANGE NOTE (bug fix): same class of issue as dominant_freq — an
% all-NaN y_profile used to silently produce a fake wavelength of
% ~1.005 BL (the first non-zero spatial-frequency bin) because MATLAB's
% max() on all-NaN input returns index 1 rather than failing. Now
% explicitly checked first.
    if all(isnan(y_profile))
        wavelength = NaN; sf = []; power = [];
        return;
    end
    N     = length(y_profile);
    ds    = s_norm(2) - s_norm(1);
    Y     = fft(y_profile - mean(y_profile, 'omitnan'));
    power = (2/N) * abs(Y(1:floor(N/2)+1)).^2;
    sf    = (0:floor(N/2)) / (N*ds);
    valid = sf > 0;
    if ~any(valid), wavelength = NaN; return; end
    power_valid = power(valid);
    sf_valid    = sf(valid);
    [~, idx]  = max(power_valid);
    f_dom     = sf_valid(idx);
    wavelength = 1 / f_dom;
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