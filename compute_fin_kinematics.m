function fin = compute_fin_kinematics(csvPath, rootName, tipName, fps, min_freq, body_speed_BL_s)
% COMPUTE_FIN_KINEMATICS  Compute 3-D kinematics of a fin from root and tip landmarks.
%
%   fin = compute_fin_kinematics(csvPath, rootName, tipName, fps)
%   fin = compute_fin_kinematics(csvPath, rootName, tipName, fps, min_freq, body_speed_BL_s)
%
%   INPUTS
%     csvPath   - path to CSV with landmark columns.
%                 Supported: named format (landmark_X/Y/Z)
%                            numbered format (pt1_X/Y/Z, pt2_X/Y/Z, ...)
%     rootName  - base name of the fin root landmark, e.g. 'Rpectbase' or 'pt1'
%                 May be a compound average string: 'Rpectbase+Lpectbase'
%     tipName   - base name of the fin tip landmark, e.g. 'Rpecttip' or 'pt2'
%     fps       - frames per second (scalar)
%     min_freq  - OPTIONAL, Hz floor for fin-beat frequency detection
%                 (default 0.5). Same purpose as compute_kinematics' min_freq —
%                 set it based on a visual beat count for this trial, not
%                 blindly left at the default (see prior evaluation notes:
%                 a too-high floor can make the FFT lock onto a harmonic
%                 instead of the true fundamental).
%     body_speed_BL_s - OPTIONAL, mean forward body speed in BL/s (e.g.
%                 ext.mean_speed_BL_s from compute_body_extended.m) — used
%                 to compute fin stride length. If omitted, stride_length
%                 fields are NaN.
%
%   OUTPUT  fin — struct with fields (existing fields unchanged, NEW ones marked):
%
%   --- Raw data ---
%     .frames        [nFrames x 1]  frame indices
%     .valid         [nFrames x 1]  logical, true where both root & tip are present
%     .root_xyz      [nFrames x 3]  absolute root position
%     .tip_xyz       [nFrames x 3]  absolute tip position
%     .vx/vy/vz      [nFrames x 1]  fin vector components (tip - root)
%
%   --- Euler angles (degrees) ---
%     .yaw / .pitch / .roll   (see original docstring for definitions —
%                              unchanged)
%
%   --- Derived kinematics ---
%     .fin_length, .mean_length, .std_length
%     .step_dist, .cum_dist, .total_dist, .tip_speed
%     .d_yaw / .d_pitch / .d_roll, .ang_vel
%
%   --- NEW: Fin beat frequency ---
%     .fin_freq_Hz          dominant oscillation frequency of the fin's
%                            yaw angle over time (the fin-beat rate,
%                            analogous to tail-beat frequency)
%     .fin_freq_pitch_Hz    same, but from the pitch signal — useful when
%                            the fin's dominant motion is more of a
%                            flap (pitch) than a sweep (yaw); compare the
%                            two and use whichever matches your visual
%                            beat count.
%
%   --- NEW: Stride length / duration ---
%     .stride_length_BL     forward body distance traveled per fin-beat
%                            cycle = body_speed_BL_s / fin_freq_Hz
%                            (NaN if body_speed_BL_s not supplied, or if
%                            fin_freq_Hz is NaN/zero)
%     .stride_duration_s    period of one fin-beat cycle = 1/fin_freq_Hz
%
%   --- Summary statistics (valid frames only) ---
%     (unchanged from original — see original docstring)
%
%   REQUIRED TOOLBOXES: none (uses only base MATLAB)

    if nargin < 5 || isempty(min_freq), min_freq = 0.5; end
    if nargin < 6, body_speed_BL_s = NaN; end

    %% ---- Load CSV ----
    opts = detectImportOptions(csvPath);
    opts.VariableNamingRule = 'preserve';
    T    = readtable(csvPath, opts);
    cols = T.Properties.VariableNames;
    nFrames = height(T);
    frames  = (1:nFrames)';

    %% ---- Extract root and tip XYZ ----
    [root_xyz, root_err] = extract_xyz(T, cols, rootName, nFrames);
    [tip_xyz,  tip_err]  = extract_xyz(T, cols, tipName,  nFrames);

    if ~isempty(root_err)
        error('compute_fin_kinematics: root point "%s" — %s', rootName, root_err);
    end
    if ~isempty(tip_err)
        error('compute_fin_kinematics: tip point "%s" — %s', tipName, tip_err);
    end

    %% ---- Fin vector (tip relative to root) ----
    vx = tip_xyz(:,1) - root_xyz(:,1);
    vy = tip_xyz(:,2) - root_xyz(:,2);
    vz = tip_xyz(:,3) - root_xyz(:,3);

    valid = ~(isnan(vx) | isnan(vy) | isnan(vz));

    if sum(valid) < 2
        error('compute_fin_kinematics: fewer than 2 valid frames for %s -> %s (%d/%d valid).', ...
              rootName, tipName, sum(valid), nFrames);
    end

    %% ---- Fin length ----
    fin_length = sqrt(vx.^2 + vy.^2 + vz.^2);

    %% ---- Euler angles ----
    yaw = atan2d(vy, vx);
    pitch = atan2d(-vz, sqrt(vx.^2 + vy.^2));
    roll = atan2d(vz, vy);

    %% ---- Tip travel ----
    step_dist = zeros(nFrames, 1);
    for f = 2:nFrames
        if valid(f) && valid(f-1)
            step_dist(f) = sqrt((tip_xyz(f,1)-tip_xyz(f-1,1))^2 + ...
                                (tip_xyz(f,2)-tip_xyz(f-1,2))^2 + ...
                                (tip_xyz(f,3)-tip_xyz(f-1,3))^2);
        end
    end
    cum_dist   = cumsum(step_dist);
    total_dist = cum_dist(end);
    tip_speed  = step_dist * fps;

    %% ---- Angular rates (deg/s) ----
    d_yaw   = [0; diff(yaw)]   * fps;
    d_pitch = [0; diff(pitch)] * fps;
    d_roll  = [0; diff(roll)]  * fps;
    ang_vel = sqrt(d_yaw.^2 + d_pitch.^2 + d_roll.^2);

    %% ---- NEW: Fin beat frequency ----
    % Uses the same NaN-safe dominant_freq logic established in
    % compute_kinematics.m: degenerate (all-NaN / zero-variance) input
    % returns NaN rather than a fabricated frequency.
    yaw_filled   = fill_nan_local(yaw);
    pitch_filled = fill_nan_local(pitch);
    fin_freq_Hz       = dominant_freq_local(yaw_filled,   fps, min_freq);
    fin_freq_pitch_Hz = dominant_freq_local(pitch_filled, fps, min_freq);

    %% ---- NEW: Stride length / duration ----
    stride_duration_s = NaN;
    stride_length_BL  = NaN;
    if ~isnan(fin_freq_Hz) && fin_freq_Hz > 0
        stride_duration_s = 1 / fin_freq_Hz;
        if ~isnan(body_speed_BL_s)
            stride_length_BL = body_speed_BL_s / fin_freq_Hz;
        end
    end

    %% ---- Summary statistics over valid frames ----
    yaw_v       = yaw(valid);
    pitch_v     = pitch(valid);
    roll_v      = roll(valid);
    % CHANGE NOTE (bug fix): frame 1 has no prior frame, so its step_dist
    % is 0 and including it in the speed stats dragged mean_speed down by
    % 1/nFrames. Only frames with a valid previous frame carry a real
    % step speed.
    has_prev    = [false; valid(1:end-1)];
    spd_v       = tip_speed(valid & has_prev);
    ang_vel_v   = ang_vel(valid);
    len_v       = fin_length(valid);

    %% ---- Pack output ----
    fin.rootName   = rootName;
    fin.tipName    = tipName;
    fin.fps        = fps;
    fin.frames     = frames;
    fin.valid      = valid;
    fin.root_xyz   = root_xyz;
    fin.tip_xyz    = tip_xyz;
    fin.vx = vx;  fin.vy = vy;  fin.vz = vz;

    fin.fin_length  = fin_length;
    fin.mean_length = mean(len_v,  'omitnan');
    fin.std_length  = std(len_v,   'omitnan');

    fin.yaw    = yaw;
    fin.pitch  = pitch;
    fin.roll   = roll;

    fin.step_dist  = step_dist;
    fin.cum_dist   = cum_dist;
    fin.total_dist = total_dist;
    fin.tip_speed  = tip_speed;

    fin.d_yaw   = d_yaw;
    fin.d_pitch = d_pitch;
    fin.d_roll  = d_roll;
    fin.ang_vel = ang_vel;

    fin.fin_freq_Hz       = fin_freq_Hz;        % NEW
    fin.fin_freq_pitch_Hz = fin_freq_pitch_Hz;  % NEW
    fin.stride_duration_s = stride_duration_s;  % NEW
    fin.stride_length_BL  = stride_length_BL;   % NEW

    fin.mean_yaw    = mean(yaw_v);     fin.std_yaw    = std(yaw_v);
    fin.range_yaw   = range(yaw_v);
    fin.mean_pitch  = mean(pitch_v);   fin.std_pitch  = std(pitch_v);
    fin.range_pitch = range(pitch_v);
    fin.mean_roll   = mean(roll_v);    fin.std_roll   = std(roll_v);
    fin.range_roll  = range(roll_v);

    fin.mean_speed   = mean(spd_v,     'omitnan');
    fin.std_speed    = std(spd_v,      'omitnan');
    fin.peak_speed   = max(spd_v);
    fin.mean_ang_vel = mean(ang_vel_v, 'omitnan');
    fin.std_ang_vel  = std(ang_vel_v,  'omitnan');
    fin.peak_ang_vel = max(ang_vel_v);

    % Tracked window: rows where the fin has any data — excludes
    % leading/trailing all-empty padding rows common in exported CSVs.
    n_fin_window = sum(any(~isnan(root_xyz) | ~isnan(tip_xyz), 2));

    fin.n_valid   = sum(valid);
    fin.n_frames  = nFrames;
    fin.n_tracked_frames = n_fin_window;
    fin.pct_valid = 100 * sum(valid) / max(n_fin_window, 1);

    fprintf('Fin: %s -> %s | %d/%d tracked frames valid | length=%.3f\xB1%.3f\n', ...
            rootName, tipName, fin.n_valid, n_fin_window, fin.mean_length, fin.std_length);
    fprintf('  Yaw  : mean=%.1f SD=%.1f range=%.1f deg\n', fin.mean_yaw,   fin.std_yaw,   fin.range_yaw);
    fprintf('  Pitch: mean=%.1f SD=%.1f range=%.1f deg\n', fin.mean_pitch, fin.std_pitch, fin.range_pitch);
    fprintf('  Roll : mean=%.1f SD=%.1f range=%.1f deg\n', fin.mean_roll,  fin.std_roll,  fin.range_roll);
    fprintf('  Tip dist=%.4f  peak speed=%.4f/s  peak ang_vel=%.2f deg/s\n', ...
            fin.total_dist, fin.peak_speed, fin.peak_ang_vel);
    fprintf('  Fin freq (yaw)=%s Hz  (pitch)=%s Hz  stride=%s BL  stride duration=%s s\n', ...
            fmt_val(fin_freq_Hz), fmt_val(fin_freq_pitch_Hz), ...
            fmt_val(stride_length_BL), fmt_val(stride_duration_s));
end


% =========================================================================
%  LOCAL HELPERS
% =========================================================================

function s = fmt_val(v)
    if isnan(v), s = 'NaN'; else, s = sprintf('%.4f', v); end
end

function y = fill_nan_local(y)
    t = (1:length(y))';
    good = ~isnan(y);
    if sum(good) < 2, y(:) = NaN; return; end
    y(~good) = interp1(t(good), y(good), t(~good), 'linear', 'extrap');
end

function f_dom = dominant_freq_local(y, fs, min_freq)
% Same NaN-safe convention as compute_kinematics.m's fixed dominant_freq:
% degenerate input -> NaN, never a fabricated bin.
    if all(isnan(y)) || (max(y) - min(y)) < eps
        f_dom = NaN; return;
    end
    N = length(y);
    Y = fft(y - mean(y));
    power = (2/N) * abs(Y(1:floor(N/2)+1)).^2;
    freqs = fs * (0:floor(N/2)) / N;
    valid = freqs >= min_freq;
    if ~any(valid), f_dom = NaN; return; end
    power_valid = power(valid);
    freqs_valid = freqs(valid);
    [~, idx] = max(power_valid);
    f_dom = freqs_valid(idx);
end

function [xyz, err] = extract_xyz(T, cols, base, nFrames)
% Extract X/Y/Z for one landmark (handles compound 'A+B' averages).
    err = '';
    parts = strtrim(strsplit(strtrim(base), '+'));

    if numel(parts) > 1
        acc = zeros(nFrames, 3);
        cnt = zeros(nFrames, 3);
        for p = 1:numel(parts)
            [sub, sub_err] = extract_xyz(T, cols, parts{p}, nFrames);
            if ~isempty(sub_err), xyz = []; err = sub_err; return; end
            mask = ~isnan(sub);
            acc(mask) = acc(mask) + sub(mask);
            cnt(mask) = cnt(mask) + 1;
        end
        cnt(cnt == 0) = NaN;
        xyz = acc ./ cnt;
        return;
    end

    xyz = NaN(nFrames, 3);
    dim_labels = {'X','Y','Z'};
    for d = 1:3
        candidates = {[base '_' dim_labels{d}], [base '_' lower(dim_labels{d})]};
        found = false;
        for v = 1:numel(candidates)
            if ismember(candidates{v}, cols)
                col_data = T.(candidates{v});
                % CHANGE NOTE (bug fix): if EVERY value in this column is
                % missing (e.g. a fully-occluded landmark that is "NA" on
                % every row), MATLAB's readtable can't infer it's numeric
                % and imports it as text (cell/char) instead of double.
                % Assigning that straight into a double array used to
                % throw "Conversion to double from cell is not possible."
                % Detect and convert explicitly instead — str2double()
                % turns unparseable text (including "NA") into NaN, which
                % is exactly the behavior we want: a fully-occluded point
                % becomes a clean all-NaN column, not a crash.
                if iscell(col_data) || isstring(col_data)
                    col_data = str2double(col_data);
                end
                xyz(:,d) = col_data;
                found = true; break;
            end
        end
        if ~found
            pat = sprintf('^%s_[%s%s]$', regexptranslate('escape',base), dim_labels{d}, lower(dim_labels{d}));
            idx = find(~cellfun(@isempty, regexpi(cols, pat)), 1);
            if ~isempty(idx)
                col_data = T.(cols{idx});
                if iscell(col_data) || isstring(col_data)
                    col_data = str2double(col_data);
                end
                xyz(:,d) = col_data;
            else
                xyz = [];
                nearby = strjoin(cols(contains(cols, base, 'IgnoreCase',true)), ', ');
                if isempty(nearby)
                    err = sprintf('Column "%s_%s" not found.', base, dim_labels{d});
                else
                    err = sprintf('Column "%s_%s" not found. Similar columns: %s', base, dim_labels{d}, nearby);
                end
                return;
            end
        end
    end

    % NEW: if the resulting column is 100% NaN, warn clearly rather than
    % letting it silently propagate into "fewer than 2 valid frames" deep
    % inside compute_fin_kinematics with no context about WHY.
    for d = 1:3
        if all(isnan(xyz(:,d)))
            warning(['extract_xyz: "%s_%s" is 100%% missing in this file (fully ' ...
                     'occluded/untracked landmark) — likely the far side of the ' ...
                     'animal from the camera. Consider using the corresponding ' ...
                     'left/near-side point instead.'], base, dim_labels{d});
        end
    end
end
