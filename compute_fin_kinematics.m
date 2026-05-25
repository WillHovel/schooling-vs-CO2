function fin = compute_fin_kinematics(csvPath, rootName, tipName, fps)
% COMPUTE_FIN_KINEMATICS  Compute 3-D kinematics of a fin from root and tip landmarks.
%
%   fin = compute_fin_kinematics(csvPath, rootName, tipName, fps)
%
%   INPUTS
%     csvPath   - path to CSV with landmark columns.
%                 Supported: named format (landmark_X/Y/Z)
%                            numbered format (pt1_X/Y/Z, pt2_X/Y/Z, ...)
%     rootName  - base name of the fin root landmark, e.g. 'Rpectbase' or 'pt1'
%                 May be a compound average string: 'Rpectbase+Lpectbase'
%     tipName   - base name of the fin tip landmark, e.g. 'Rpecttip' or 'pt2'
%     fps       - frames per second (scalar)
%
%   OUTPUT  fin — struct with fields:
%
%   --- Raw data ---
%     .frames        [nFrames x 1]  frame indices
%     .valid         [nFrames x 1]  logical, true where both root & tip are present
%     .root_xyz      [nFrames x 3]  absolute root position
%     .tip_xyz       [nFrames x 3]  absolute tip position
%     .vx/vy/vz      [nFrames x 1]  fin vector components (tip - root)
%
%   --- Euler angles (degrees) ---
%     .yaw           azimuth of fin vector in XY plane  (atan2(vy, vx))
%     .pitch         elevation out of XY plane          (atan2(-vz, sqrt(vx^2+vy^2)))
%     .roll          rotation about the fin's long axis, estimated as the angle
%                    between the fin vector projected into the YZ plane and the
%                    world +Z axis — i.e. how much the fin tilts dorsoventrally
%                    relative to horizontal.  Range: -180 to +180 deg.
%                    NOTE: true roll (rotation about the fin's own axis) is
%                    geometrically undefined from a single vector.  This metric
%                    is the closest physically meaningful proxy: it captures
%                    fin dorsoventral cant about the body long axis.
%
%   --- Derived kinematics ---
%     .fin_length     [nFrames x 1]  Euclidean root-to-tip distance per frame
%     .mean_length / .std_length     mean/SD fin length across valid frames
%     .step_dist      [nFrames x 1]  tip displacement frame-to-frame (mm equiv)
%     .cum_dist       [nFrames x 1]  cumulative tip distance
%     .total_dist     scalar         total tip path length
%     .tip_speed      [nFrames x 1]  tip speed (units/s)
%     .d_yaw / d_pitch / d_roll      [nFrames x 1] angular rates (deg/s)
%     .ang_vel        [nFrames x 1]  combined angular velocity magnitude (deg/s)
%
%   --- Summary statistics (valid frames only) ---
%     .mean_yaw   / .std_yaw   / .range_yaw
%     .mean_pitch / .std_pitch / .range_pitch
%     .mean_roll  / .std_roll  / .range_roll
%     .mean_speed / .std_speed / .peak_speed
%     .mean_ang_vel / .std_ang_vel / .peak_ang_vel
%     .n_valid / .n_frames / .pct_valid
%
%   REQUIRED TOOLBOXES: none (uses only base MATLAB)

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
    % YAW: azimuth in the XY (horizontal) plane — how far forward/backward the
    %      fin sweeps.  0 deg = pointing in +X; 90 deg = pointing in +Y.
    yaw = atan2d(vy, vx);

    % PITCH: elevation above/below the XY plane — how much the fin is raised or
    %        depressed.  +90 = pointing straight up (+Z), -90 = straight down.
    pitch = atan2d(-vz, sqrt(vx.^2 + vy.^2));

    % ROLL: dorsoventral cant of the fin about the body's long axis (X-axis).
    %       Defined as atan2(vz, vy) — the angle of the fin vector in the YZ
    %       plane (the cross-section plane of the body).
    %       0 deg = fin pointing in +Y (lateral); +90 = fin pointing up (+Z).
    %       This is the closest proxy to "roll" achievable from a single vector.
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

    %% ---- Summary statistics over valid frames ----
    yaw_v       = yaw(valid);
    pitch_v     = pitch(valid);
    roll_v      = roll(valid);
    spd_v       = tip_speed(valid);
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

    fin.n_valid  = sum(valid);
    fin.n_frames = nFrames;
    fin.pct_valid = 100 * sum(valid) / nFrames;

    fprintf('Fin: %s -> %s | %d/%d valid frames | length=%.3f±%.3f\n', ...
            rootName, tipName, fin.n_valid, nFrames, fin.mean_length, fin.std_length);
    fprintf('  Yaw  : mean=%.1f SD=%.1f range=%.1f deg\n', fin.mean_yaw,   fin.std_yaw,   fin.range_yaw);
    fprintf('  Pitch: mean=%.1f SD=%.1f range=%.1f deg\n', fin.mean_pitch, fin.std_pitch, fin.range_pitch);
    fprintf('  Roll : mean=%.1f SD=%.1f range=%.1f deg\n', fin.mean_roll,  fin.std_roll,  fin.range_roll);
    fprintf('  Tip dist=%.4f  peak speed=%.4f/s  peak ang_vel=%.2f deg/s\n', ...
            fin.total_dist, fin.peak_speed, fin.peak_ang_vel);
end


% =========================================================================
%  LOCAL HELPERS
% =========================================================================

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
                xyz(:,d) = T.(candidates{v});
                found = true; break;
            end
        end
        if ~found
            pat = sprintf('^%s_[%s%s]$', regexptranslate('escape',base), dim_labels{d}, lower(dim_labels{d}));
            idx = find(~cellfun(@isempty, regexpi(cols, pat)), 1);
            if ~isempty(idx)
                xyz(:,d) = T.(cols{idx});
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
end