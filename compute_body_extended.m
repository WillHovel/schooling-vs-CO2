function ext = compute_body_extended(fish_points, fps, kine, roll_pair, flow_BL_s)
% COMPUTE_BODY_EXTENDED  Additional whole-body kinematics not covered by
%                        compute_kinematics.m: body angle, speed, stride
%                        length, head elevation, and roll (if available).
%
%   ext = compute_body_extended(fish_points, fps, kine, roll_pair)
%   ext = compute_body_extended(fish_points, fps, kine, roll_pair, flow_BL_s)
%
%   INPUTS
%     fish_points - struct from transform_fish() — must have .points
%                   (raw, untransformed coords), .X/.Y[/.Z], and
%                   .transform_params (added by the updated transform_fish).
%     fps         - frames per second.
%     kine        - output of compute_kinematics(fish_points, fps, min_freq)
%                   for the SAME fish — used to get tail_TBF for stride
%                   length and Strouhal. Pass [] to skip stride length.
%     roll_pair   - OPTIONAL 1x2 cell of point_names, e.g.
%                   {'LPectBase','RPectBase'}, giving a left/right pair
%                   used to estimate roll. Roll cannot be computed from a
%                   single midline — it needs two laterally-symmetric
%                   points to tell which way "up" tilts. If omitted or
%                   the named points aren't found, .roll_deg is NaN and
%                   .roll_available is false.
%     flow_BL_s   - OPTIONAL tank/flume flow speed in BL/s, SIGNED:
%                     + = flow opposes the fish (upstream station-holding;
%                         through-water speed = ground speed + flow)
%                     - = flow assists the fish (downstream;
%                         through-water speed = |ground speed - |flow||)
%                    0 or omitted = no flow correction (ground speed only).
%                   This must be converted to BL/s by the CALLER (flow in
%                   cm/s divided by body length in cm, etc.) — this
%                   function works entirely in BL units.
%
%   OUTPUT  ext — struct (one per animal) with fields:
%
%   --- Body angle (heading in the raw/world coordinate frame) ---
%     .body_angle_deg      [nFrames x 1]  per-frame swim-axis angle (deg),
%                           i.e. how the fish's own heading changes over
%                           time relative to the camera/world frame — NOT
%                           the same as lateral undulation; this is the
%                           slope of the fitted body axis before rotation.
%                           Continuous mod-180-unwrapped heading series
%                           (mean/std relative to the first valid frame).
%                           NaN for pre-transformed (CURVES) data, whose
%                           heading information no longer exists.
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
%   --- Through-water speed & Strouhal (needs flow_BL_s input) ---
%     .flow_BL_s             the signed flow speed used (0 = none)
%     .flow_orientation      'against' / 'with' / 'none'
%     .speed_through_water_BL_s   [nFrames x 1] ground speed corrected for
%                                  flow (see flow_BL_s sign convention above)
%     .mean_speed_through_water_BL_s / .std_... / .peak_...
%     .tail_amp_pp_BL        peak-to-peak trailing-edge (tail point) lateral
%                             excursion in BL = 2*sqrt(2)*std(tail Y) — the
%                             standard sinusoid amplitude estimator
%     .strouhal               St = tail_TBF * tail_amp_pp_BL / U, where
%                             U = mean through-water speed if flow is set,
%                             otherwise mean ground speed. Dimensionless.
%                             (Triantafyllou convention: peak-to-peak
%                             trailing-edge amplitude.) NaN when tail_TBF,
%                             amplitude, or U is missing/zero.
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

    if nargin < 5 || isempty(flow_BL_s), flow_BL_s = 0; end
    if ~isscalar(flow_BL_s), flow_BL_s = flow_BL_s(1); end

    nFish = numel(fish_points);
    ext(nFish) = struct();

    for fi = 1:nFish
        pts     = fish_points(fi).points;     % [nFrames x nPoints x nDims]
        nFrames = size(pts, 1);
        nPoints = size(pts, 2);
        has_z   = isfield(fish_points(fi),'has_z') && fish_points(fi).has_z;
        middle_idx = 2:nPoints-1;
        % Same degenerate-fit guard as transform_fish: with a single
        % middle point, polyfit returns the meaningless slope y/x — fall
        % back to the head-to-tail chord (endpoints) instead.
        if numel(middle_idx) >= 2
            fit_idx = middle_idx;
        else
            fit_idx = [1, nPoints];
        end

        % ================================================================
        % 1. BODY ANGLE — recompute the same middle-point line fit
        %    transform_fish uses, but report the angle itself instead of
        %    using it to rotate. This is the fish's heading in the raw
        %    (world/camera) coordinate frame.
        %
        %    CHANGE NOTE (bug fix): atand() returns a line ORIENTATION in
        %    (-90, 90] — angles differing by 180 deg are the same line. A
        %    fish that slowly turns past +/-90 deg used to produce a fake
        %    ~180 deg jump in the per-frame trace (and garbage mean/std/
        %    range/angular-velocity). Now the trace is unwrapped mod 180
        %    deg, i.e. each step is forced into (-90, 90], so the output
        %    is a continuous heading series (mean/std are relative to the
        %    first valid frame's value).
        %
        %    For pre-transformed data (CURVES) the heading information no
        %    longer exists (points are already in the body frame), so the
        %    angle is left NaN rather than reporting a meaningless ~0.
        % ================================================================
        pre_xformed = isfield(fish_points(fi),'pre_transformed') && fish_points(fi).pre_transformed;

        body_angle_deg = NaN(nFrames, 1);
        if ~pre_xformed
            for f = 1:nFrames
                x_mid = squeeze(pts(f, fit_idx, 1));
                y_mid = squeeze(pts(f, fit_idx, 2));
                if any(isnan(x_mid)) || any(isnan(y_mid)), continue; end
                b = polyfit(x_mid, y_mid, 1);
                body_angle_deg(f) = atand(b(1));   % slope -> degrees, range (-90, 90]
            end
        end

        % Mod-180 unwrap: force each step into (-90, 90] so the heading
        % trace is continuous across the atand() wraparound.
        body_angle_unwrapped = NaN(nFrames, 1);
        last_ang = NaN;
        for f = 1:nFrames
            a = body_angle_deg(f);
            if isnan(a), continue; end
            if isnan(last_ang)
                body_angle_unwrapped(f) = a;
            else
                d = mod(a - last_ang + 90, 180) - 90;
                body_angle_unwrapped(f) = last_ang + d;
            end
            last_ang = body_angle_unwrapped(f);
        end

        angular_velocity_deg_s = [NaN; diff(body_angle_unwrapped)] * fps;

        valid_ang = ~isnan(body_angle_unwrapped);
        if any(valid_ang)
            mean_body_angle  = mean(body_angle_unwrapped(valid_ang));
            std_body_angle   = std(body_angle_unwrapped(valid_ang));
            range_body_angle = range(body_angle_unwrapped(valid_ang));
        else
            mean_body_angle = NaN; std_body_angle = NaN; range_body_angle = NaN;
        end

        % ================================================================
        % 2. SPEED — frame-to-frame displacement of the body centroid, in
        %    raw units, converted to BL/s using this frame's raw body
        %    length (bl_per_frame, from transform_fish).
        %
        %    CHANGE NOTE (bug fix): pre-transformed data (CURVES, Format E)
        %    never gets a .bl_per_frame field because transform_fish skips
        %    it — this used to throw "Unrecognized field name" and take
        %    down the whole run. Such data is ALREADY in BL units (X from
        %    0=head to ~1=tail), so the per-frame body length is exactly 1
        %    and speed is simply centroid displacement x fps.
        % ================================================================
        centroid = squeeze(mean(pts(:,:,1:2), 2, 'omitnan'));   % [nFrames x 2]

        if isfield(fish_points(fi),'pre_transformed') && fish_points(fi).pre_transformed
            bl_pf = ones(nFrames, 1);    % CURVES data already in BL units
        elseif isfield(fish_points(fi), 'bl_per_frame') && ~isempty(fish_points(fi).bl_per_frame)
            bl_pf = fish_points(fi).bl_per_frame;   % [nFrames x 1], raw units
        else
            warning(['compute_body_extended: %s has no bl_per_frame — speed/stride ' ...
                     'will be NaN. Run transform_fish on this animal first.'], ...
                     fish_points(fi).name);
            bl_pf = NaN(nFrames, 1);
        end

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
        % 2b. THROUGH-WATER SPEED — ground speed corrected by a known
        %     tank/flume flow. Sign convention (documented above):
        %       flow > 0  (against fish):  U_tw = U_ground + flow
        %       flow < 0  (with fish):     U_tw = |U_ground - |flow||
        %     ground speed here is a MAGNITUDE (centroid distance per
        %     second), so this assumes the fish's path is aligned with the
        %     flow axis — the standard flume station-holding assumption.
        % ================================================================
        if isfinite(flow_BL_s) && flow_BL_s > 0
            speed_tw_BL_s = speed_BL_s + flow_BL_s;
            flow_orientation = 'against';
        elseif isfinite(flow_BL_s) && flow_BL_s < 0
            speed_tw_BL_s = abs(speed_BL_s + flow_BL_s);
            flow_orientation = 'with';
        else
            speed_tw_BL_s = speed_BL_s;
            flow_orientation = 'none';
        end

        valid_tw = ~isnan(speed_tw_BL_s);
        if any(valid_tw)
            mean_tw = mean(speed_tw_BL_s(valid_tw));
            std_tw  = std(speed_tw_BL_s(valid_tw));
            peak_tw = max(speed_tw_BL_s(valid_tw));
        else
            mean_tw = NaN; std_tw = NaN; peak_tw = NaN;
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
        % 3b. STROUHAL NUMBER — St = tail_TBF * A_pp / U
        %     A_pp = peak-to-peak trailing-edge (tail point) lateral
        %     excursion in BL. For a sinusoid, amplitude A = sqrt(2)*std,
        %     so A_pp = 2*sqrt(2)*std(tail Y). U = mean through-water
        %     speed when flow is set (the physically correct reference
        %     velocity in a flume), otherwise mean ground speed.
        % ================================================================
        tail_amp_pp_BL = NaN;
        if isfield(fish_points(fi),'Y') && ~isempty(fish_points(fi).Y)
            yt = fish_points(fi).Y(:, end);
            if ~pre_xformed
                % The per-frame body-axis line fit oscillates at the beat
                % frequency because the traveling wave biases it, which
                % inflates the tail's apparent excursion (~50% on clean
                % synthetic waves). Fix: measure the tail's PERPENDICULAR
                % distance from a SMOOTHED version of the fitted axis
                % (line orientation + middle-point centroid averaged over
                % ~1.5 beat periods). The orientation is smoothed as the
                % complex vector exp(i*2*alpha) — invariant to the line's
                % +/-180 deg ambiguity — so it also survives the large
                % heading swings / mirror flips of turning or walking
                % fish. Validated: synthetic raw A_pp=0.153 vs true 0.100
                % -> 0.104; robust to slow drift and to walking-shark
                % heading swings of +/-69 deg.
                tp = fish_points(fi).transform_params;
                if numel(tp) == nFrames && all(isfield(tp, {'theta','bl'}))
                    bl_pf = [tp.bl]';
                    alpha = 2*pi - [tp.theta]';   % atan(b): raw fit-line angle
                    z = exp(2i * alpha);          % direction, +/-180-invariant
                    tbf = NaN;
                    if ~isempty(kine) && isfield(kine(fi),'tail_TBF')
                        tbf = kine(fi).tail_TBF;
                    end
                    if isfinite(tbf) && tbf > 0
                        win = round(1.5 * fps / tbf);
                    else
                        win = 2 * fps;   % no beat detected — 2 s window
                    end
                    win = max(win, 5);
                    win = min(win, nFrames);
                    if mod(win,2) == 0, win = win + 1; end
                    z_s   = movmean(z, win, 'omitnan');
                    alpha_s = 0.5 * angle(z_s);
                    x_mid_all = squeeze(mean(pts(:, fit_idx, 1), 2, 'omitnan'));
                    y_mid_all = squeeze(mean(pts(:, fit_idx, 2), 2, 'omitnan'));
                    xb_s = movmean(x_mid_all, win, 'omitnan');
                    yb_s = movmean(y_mid_all, win, 'omitnan');
                    for f = 1:nFrames
                        if isnan(yt(f)) || ~(bl_pf(f) > 0), continue; end
                        x_t = pts(f, nPoints, 1);
                        y_t = pts(f, nPoints, 2);
                        if isnan(x_t) || isnan(y_t) || isnan(alpha_s(f)) ...
                                || isnan(xb_s(f)) || isnan(yb_s(f)), continue; end
                        % perpendicular distance from the smoothed axis, in BL
                        yt(f) = ((y_t - yb_s(f))*cos(alpha_s(f)) - (x_t - xb_s(f))*sin(alpha_s(f))) / bl_pf(f);
                    end
                end
            end
            good = isfinite(yt);
            if sum(good) >= 3
                tg = (1:nFrames)';
                yt(good) = detrend(yt(good), 'linear', 'SamplePoints', tg(good));
            end
            yt = yt(isfinite(yt));
            if numel(yt) >= 3
                tail_amp_pp_BL = 2 * sqrt(2) * std(yt);
            end
        end

        strouhal = NaN;
        if ~isempty(kine) && isfield(kine(fi), 'tail_TBF') && ~isnan(kine(fi).tail_TBF) ...
                && kine(fi).tail_TBF > 0 && ~isnan(tail_amp_pp_BL) && tail_amp_pp_BL > 0 ...
                && isfinite(mean_tw) && mean_tw > 0
            strouhal = kine(fi).tail_TBF * tail_amp_pp_BL / mean_tw;
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
                range_head_pitch = range(head_pitch_deg(valid_pitch));
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
                    range_roll = range(roll_deg(valid_roll));
                end
            else
                warning(['compute_body_extended: roll_pair points "%s"/"%s" not found ' ...
                         'in %s''s point_names — roll not computed.'], ...
                         roll_pair{1}, roll_pair{2}, fish_points(fi).name);
            end
        end

        % ---- Pack ----
        ext(fi).name                    = fish_points(fi).name;
        ext(fi).body_angle_deg          = body_angle_unwrapped;
        ext(fi).mean_body_angle_deg     = mean_body_angle;
        ext(fi).std_body_angle_deg      = std_body_angle;
        ext(fi).range_body_angle_deg    = range_body_angle;
        ext(fi).angular_velocity_deg_s  = angular_velocity_deg_s;

        ext(fi).speed_BL_s      = speed_BL_s;
        ext(fi).mean_speed_BL_s = mean_speed;
        ext(fi).std_speed_BL_s  = std_speed;
        ext(fi).peak_speed_BL_s = peak_speed;

        ext(fi).stride_length_BL = stride_length_BL;

        ext(fi).flow_BL_s         = flow_BL_s;
        ext(fi).flow_orientation  = flow_orientation;
        ext(fi).speed_through_water_BL_s      = speed_tw_BL_s;
        ext(fi).mean_speed_through_water_BL_s = mean_tw;
        ext(fi).std_speed_through_water_BL_s  = std_tw;
        ext(fi).peak_speed_through_water_BL_s = peak_tw;
        ext(fi).tail_amp_pp_BL    = tail_amp_pp_BL;
        ext(fi).strouhal          = strouhal;

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
        if isfinite(flow_BL_s) && flow_BL_s ~= 0
            fprintf('  flow=%.3f BL/s (%s) -> through-water speed=%.3f\xB1%.3f BL/s\n', ...
                    flow_BL_s, flow_orientation, mean_tw, std_tw);
        end
        if isfinite(strouhal)
            fprintf('  tail amp (p-p)=%.4f BL  Strouhal=%.3f (U=%s)\n', ...
                    tail_amp_pp_BL, strouhal, flow_orientation);
        end
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
