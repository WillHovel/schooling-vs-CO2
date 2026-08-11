function ext = compute_body_extended(fish_points, fps, kine, roll_pair)
% COMPUTE_BODY_EXTENDED  Additional whole-body kinematics not covered by
%                        compute_kinematics.m: body angle, speed, stride
%                        length, head elevation, and roll (if available).
%
%   ext = compute_body_extended(fish_points, fps, kine, roll_pair)
%
%   INPUTS
%     fish_points - struct from transform_fish() — must have .points
%                   (raw, untransformed coords), .X/.Y[/.Z], and
%                   .transform_params (added by the updated transform_fish).
%     fps         - frames per second.
%     kine        - output of compute_kinematics(fish_points, fps, min_freq)
%                   for the SAME fish — used to get tail_TBF for stride
%                   length. Pass [] to skip stride length.
%     roll_pair   - OPTIONAL 1x2 cell of point_names, e.g.
%                   {'LPectBase','RPectBase'}, giving a left/right pair
%                   used to estimate roll. Roll cannot be computed from a
%                   single midline — it needs two laterally-symmetric
%                   points to tell which way "up" tilts. If omitted or
%                   the named points aren't found, .roll_deg is NaN and
%                   .roll_available is false.
%
%   OUTPUT  ext — struct (one per animal) with fields:
%
%   --- Body angle (heading in the raw/world coordinate frame) ---
%     .body_angle_deg      [nFrames x 1]  per-frame swim-axis angle (deg),
%                           i.e. how the fish's own heading changes over
%                           time relative to the camera/world frame — NOT
%                           the same as lateral undulation; this is the
%                           slope of the fitted body axis before rotation.
%     .mean_body_angle_deg / .std_body_angle_deg / .range_body_angle_deg
%     .angular_velocity_deg_s   [nFrames x 1]  turning rate (deg/s)
%
%   --- Speed ---
%     .speed_BL_s           [nFrames x 1]  forward speed of the body
%                            centroid, in body-lengths/second (BL/s is the
%                            standard unit in swimming kinematics — avoids
%                            needing absolute camera calibration).
%     .mean_speed_BL_s / .std_speed_BL_s / .peak_speed_BL_s
%
%   --- Stride length (needs kine for tail_TBF) ---
%     .stride_length_BL      mean forward distance traveled per tail-beat
%                             cycle = mean_speed_BL_s / tail_TBF (BL/cycle)
%
%   --- Head elevation ---
%     .head_pitch_deg        [nFrames x 1]  angle of the head-to-next-point
%                             vector relative to horizontal (deg) — this is
%                             what most kinematics papers mean by "head
%                             elevation angle." Only computed if 3-D data
%                             (has_z) is available.
%     .mean_head_pitch_deg / .std_head_pitch_deg / .range_head_pitch_deg
%     .head_Z_raw             [nFrames x 1]  raw vertical (Z) position of
%                              the head point, relative to its own first
%                              valid frame — a simpler "is it swimming up
%                              or down in the water column" trace.
%
%   --- Roll (only if roll_pair given and both points found) ---
%     .roll_available         logical
%     .roll_deg                [nFrames x 1]  or NaN(nFrames,1) if unavailable
%     .mean_roll_deg / .std_roll_deg / .range_roll_deg

    nFish = numel(fish_points);
    ext(nFish) = struct();

    for fi = 1:nFish
        pts     = fish_points(fi).points;     % [nFrames x nPoints x nDims]
        nFrames = size(pts, 1);
        nPoints = size(pts, 2);
        has_z   = isfield(fish_points(fi),'has_z') && fish_points(fi).has_z;
        middle_idx = 2:nPoints-1;

        % ================================================================
        % 1. BODY ANGLE — recompute the same middle-point line fit
        %    transform_fish uses, but report the angle itself instead of
        %    using it to rotate. This is the fish's heading in the raw
        %    (world/camera) coordinate frame.
        % ================================================================
        body_angle_deg = NaN(nFrames, 1);
        for f = 1:nFrames
            x_mid = squeeze(pts(f, middle_idx, 1));
            y_mid = squeeze(pts(f, middle_idx, 2));
            if any(isnan(x_mid)) || any(isnan(y_mid)), continue; end
            b = polyfit(x_mid, y_mid, 1);
            body_angle_deg(f) = atand(b(1));   % slope -> degrees, range (-90, 90]
        end

        % Unwrap-free angular velocity (deg/s) — fine for a (-90,90] range
        % since real frame-to-frame turning is small relative to 90 deg.
        angular_velocity_deg_s = [NaN; diff(body_angle_deg)] * fps;

        valid_ang = ~isnan(body_angle_deg);
        if any(valid_ang)
            mean_body_angle = mean(body_angle_deg(valid_ang));
            std_body_angle  = std(body_angle_deg(valid_ang));
            range_body_angle = max(body_angle_deg(valid_ang)) - min(body_angle_deg(valid_ang));
        else
            mean_body_angle = NaN; std_body_angle = NaN; range_body_angle = NaN;
        end

        % ================================================================
        % 2. SPEED — frame-to-frame displacement of the body centroid, in
        %    raw units, converted to BL/s using this frame's raw body
        %    length (bl_per_frame, from transform_fish).
        % ================================================================
        centroid = squeeze(mean(pts(:,:,1:2), 2, 'omitnan'));   % [nFrames x 2]
        bl_pf    = fish_points(fi).bl_per_frame;                 % [nFrames x 1], raw units

        speed_BL_s = NaN(nFrames, 1);
        for f = 2:nFrames
            if any(isnan(centroid(f,:))) || any(isnan(centroid(f-1,:))), continue; end
            bl_here = bl_pf(f);
            if isnan(bl_here) || bl_here <= 0
                bl_here = bl_pf(f-1);   % fall back to adjacent frame's BL if this one's missing
            end
            if isnan(bl_here) || bl_here <= 0, continue; end
            step_raw = sqrt(sum((centroid(f,:) - centroid(f-1,:)).^2));
            speed_BL_s(f) = (step_raw / bl_here) * fps;
        end

        valid_spd = ~isnan(speed_BL_s);
        if any(valid_spd)
            mean_speed = mean(speed_BL_s(valid_spd));
            std_speed  = std(speed_BL_s(valid_spd));
            peak_speed = max(speed_BL_s(valid_spd));
        else
            mean_speed = NaN; std_speed = NaN; peak_speed = NaN;
        end

        % ================================================================
        % 3. STRIDE LENGTH — mean forward distance traveled per tail-beat
        %    cycle. Needs tail_TBF from compute_kinematics.
        % ================================================================
        stride_length_BL = NaN;
        if ~isempty(kine) && isfield(kine(fi), 'tail_TBF') && ~isnan(kine(fi).tail_TBF) ...
                && kine(fi).tail_TBF > 0 && ~isnan(mean_speed)
            stride_length_BL = mean_speed / kine(fi).tail_TBF;
        end

        % ================================================================
        % 4. HEAD ELEVATION
        % ================================================================
        head_pitch_deg = NaN(nFrames, 1);
        head_Z_raw     = NaN(nFrames, 1);
        mean_head_pitch = NaN; std_head_pitch = NaN; range_head_pitch = NaN;

        if has_z && nPoints >= 2
            z_head_first_valid = NaN;
            for f = 1:nFrames
                x1 = pts(f,1,1); y1 = pts(f,1,2); z1 = pts(f,1,3);
                x2 = pts(f,2,1); y2 = pts(f,2,2); z2 = pts(f,2,3);
                if any(isnan([x1 y1 z1 x2 y2 z2])), continue; end

                horiz_dist = sqrt((x2-x1)^2 + (y2-y1)^2);
                if horiz_dist > eps
                    head_pitch_deg(f) = atand((z1 - z2) / horiz_dist);  % + = head raised
                end

                if isnan(z_head_first_valid), z_head_first_valid = z1; end
                head_Z_raw(f) = z1 - z_head_first_valid;
            end
            valid_pitch = ~isnan(head_pitch_deg);
            if any(valid_pitch)
                mean_head_pitch  = mean(head_pitch_deg(valid_pitch));
                std_head_pitch   = std(head_pitch_deg(valid_pitch));
                range_head_pitch = max(head_pitch_deg(valid_pitch)) - min(head_pitch_deg(valid_pitch));
            end
        end

        % ================================================================
        % 5. ROLL — only if a valid left/right point pair is supplied.
        % ================================================================
        roll_available = false;
        roll_deg = NaN(nFrames, 1);
        mean_roll = NaN; std_roll = NaN; range_roll = NaN;

        if nargin >= 4 && ~isempty(roll_pair) && has_z
            pn = fish_points(fi).point_names;
            iL = find(strcmpi(pn, roll_pair{1}), 1);
            iR = find(strcmpi(pn, roll_pair{2}), 1);
            if ~isempty(iL) && ~isempty(iR)
                roll_available = true;
                for f = 1:nFrames
                    yL = pts(f,iL,2); zL = pts(f,iL,3);
                    yR = pts(f,iR,2); zR = pts(f,iR,3);
                    if any(isnan([yL zL yR zR])), continue; end
                    roll_deg(f) = atan2d(zR - zL, yR - yL);
                end
                valid_roll = ~isnan(roll_deg);
                if any(valid_roll)
                    mean_roll  = mean(roll_deg(valid_roll));
                    std_roll   = std(roll_deg(valid_roll));
                    range_roll = max(roll_deg(valid_roll)) - min(roll_deg(valid_roll));
                end
            else
                warning(['compute_body_extended: roll_pair points "%s"/"%s" not found ' ...
                         'in %s''s point_names — roll not computed.'], ...
                         roll_pair{1}, roll_pair{2}, fish_points(fi).name);
            end
        end

        % ---- Pack ----
        ext(fi).name                    = fish_points(fi).name;
        ext(fi).body_angle_deg          = body_angle_deg;
        ext(fi).mean_body_angle_deg     = mean_body_angle;
        ext(fi).std_body_angle_deg      = std_body_angle;
        ext(fi).range_body_angle_deg    = range_body_angle;
        ext(fi).angular_velocity_deg_s  = angular_velocity_deg_s;

        ext(fi).speed_BL_s      = speed_BL_s;
        ext(fi).mean_speed_BL_s = mean_speed;
        ext(fi).std_speed_BL_s  = std_speed;
        ext(fi).peak_speed_BL_s = peak_speed;

        ext(fi).stride_length_BL = stride_length_BL;

        ext(fi).head_pitch_deg       = head_pitch_deg;
        ext(fi).mean_head_pitch_deg  = mean_head_pitch;
        ext(fi).std_head_pitch_deg   = std_head_pitch;
        ext(fi).range_head_pitch_deg = range_head_pitch;
        ext(fi).head_Z_raw           = head_Z_raw;

        ext(fi).roll_available  = roll_available;
        ext(fi).roll_deg        = roll_deg;
        ext(fi).mean_roll_deg   = mean_roll;
        ext(fi).std_roll_deg    = std_roll;
        ext(fi).range_roll_deg  = range_roll;

        fprintf(['%s | body_angle=%.1f\xB1%.1fdeg  speed=%.3f\xB1%.3f BL/s  ' ...
                 'stride=%s BL  head_pitch=%.1fdeg  roll=%s\n'], ...
                fish_points(fi).name, mean_body_angle, std_body_angle, ...
                mean_speed, std_speed, fmt_stride(stride_length_BL), ...
                mean_head_pitch, fmt_roll(roll_available, mean_roll));
    end
end

function s = fmt_stride(v)
    if isnan(v), s = 'NaN'; else, s = sprintf('%.4f', v); end
end

function s = fmt_roll(avail, v)
    if ~avail, s = 'N/A (no paired L/R points given)';
    elseif isnan(v), s = 'NaN';
    else, s = sprintf('%.1fdeg', v);
    end
end
