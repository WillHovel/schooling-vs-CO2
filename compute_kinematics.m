function kine = compute_kinematics(fish_points, fps, min_freq)
% COMPUTE_KINEMATICS  FFT-based kinematic analysis on transformed fish midlines.
%
%   kine = compute_kinematics(fish_points, fps, min_freq)
%
%   INPUTS
%     fish_points - struct array from load_fish_points[_named]() -> transform_fish().
%                   Must have .X and .Y [nFrames x nPoints].
%                   If .Z is present and .has_z is true, 3-D amplitude and
%                   curvature are computed in addition to the XY quantities.
%     fps         - frames per second (scalar, or [nFish x 1] vector).
%     min_freq    - minimum plausible beat frequency in Hz (e.g. 0.5).
%
%   OUTPUT  kine  - struct array (one per animal) with fields:
%
%   --- Interpolated midlines (200 points along body) ---
%     .X_interp      [nFrames x 200]
%     .Y_interp      [nFrames x 200]
%     .Z_interp      [nFrames x 200]   (3-D only)
%     .s_norm        [1 x 200]         0 = head, 1 = tail
%
%   --- Lateral (Y) amplitude ---
%     .amp_mean / .amp_std   [1 x 200]  mean/std |Y| per body position
%     .headAmp               head-region mean half-amplitude
%     .tailAmp               tail-region mean half-amplitude
%     .headTailAmpRatio
%     .minAmp / .minAmpLoc
%     .maxAmp / .maxAmpLoc
%
%   --- Dorso-ventral (Z) amplitude (3-D only) ---
%     .ampZ_mean / .ampZ_std  [1 x 200]
%     .headAmpZ / .tailAmpZ
%     .minAmpZ / .minAmpZLoc / .maxAmpZ / .maxAmpZLoc
%
%   --- Beat frequencies ---
%     .head_TBF / .tail_TBF          (Hz, from FFT of first/last point Y over time)
%     .head_fft_freq / .head_fft_power
%     .tail_fft_freq / .tail_fft_power
%     .headZ_TBF / .tailZ_TBF        (3-D only, from Z time-series)
%
%   --- Propulsive wave ---
%     .wavelength            dominant spatial wavelength (BL) from mean |Y| profile
%     .wave_spatial_freq / .wave_power
%
%   --- Curvature ---
%     .curv_mean / .curv_std  [1 x 200]  3-point geometric curvature in XY
%     .curv3d_mean / .curv3d_std         (3-D only — curvature in 3-D space)
%     .maxCurv / .maxCurvLoc
%     .maxCurv3D / .maxCurv3DLoc        (3-D only)

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

        fprintf('%s | head TBF=%.2fHz tail TBF=%.2fHz wavelength=%.3fBL maxCurv=%.3f\n', ...
                fish_points(fi).name, head_TBF, tail_TBF, wavelength, maxCurv);
    end
end


% =========================================================================
%  LOCAL HELPERS
% =========================================================================

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

    headAmp = mean(amp_mean(s_norm <= 0.05));
    tailAmp = mean(amp_mean(s_norm >= 0.95));

    [minAmp, mi] = min(amp_mean);  minLoc = s_norm(mi);
    [maxAmp, ma] = max(amp_mean);  maxLoc = s_norm(ma);
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
    [maxCurv, ci] = max(curv_mean);
    maxCurvLoc    = s_norm(ci);
end


function [f_dom, freqs, power] = dominant_freq(y, fs, min_freq)
    N     = length(y);
    Y     = fft(y - mean(y));
    power = (2/N) * abs(Y(1:floor(N/2)+1)).^2;
    freqs = fs * (0:floor(N/2)) / N;
    valid = freqs >= min_freq;
    if ~any(valid), f_dom = NaN; return; end
    [~, idx]  = max(power(valid));
    f_dom     = freqs(find(valid, idx, 'first'));
    f_dom     = freqs(find(valid));
    f_dom     = f_dom(idx);
end


function [wavelength, sf, power] = spatial_wavelength(y_profile, s_norm)
    N     = length(y_profile);
    ds    = s_norm(2) - s_norm(1);
    Y     = fft(y_profile - mean(y_profile));
    power = (2/N) * abs(Y(1:floor(N/2)+1)).^2;
    sf    = (0:floor(N/2)) / (N*ds);
    valid = sf > 0;
    if ~any(valid), wavelength = NaN; return; end
    [~, idx]  = max(power(valid));
    f_dom     = sf(find(valid));
    f_dom     = f_dom(idx);
    wavelength = 1 / f_dom;
end


function y = fill_nan(y)
    t    = (1:length(y))';
    good = ~isnan(y);
    if sum(good) < 2, y(:) = 0; return; end
    y(~good) = interp1(t(good), y(good), t(~good), 'linear', 'extrap');
end