function fish_points = transform_fish(fish_points, bl_override)
% TRANSFORM_FISH  Rotate and translate each fish's midline so the mean
%                 swimming axis is parallel to the X-axis, with the head
%                 at X = 0.
%
%   fish_points = transform_fish(fish_points)
%   fish_points = transform_fish(fish_points, bl_override)
%
%   bl_override — OPTIONAL known body length (same units as the CSV
%                 coordinates), scalar (applied to all fish) or 1xN per
%                 fish; empty/NaN = auto. Overrides the per-frame
%                 head-to-tail distance used to normalize .X/.Y[/.Z] and
%                 to fill .bl_per_frame. Use this when the first/last
%                 tracked points are not the true head/tail (e.g. Format C
%                 numbered files whose endpoints are fin bases) or when
%                 you have a calibrated body length.
%
%   Works with any number of points (>= 3) and with 2-D or 3-D data.
%   The struct schema is identical to the output of load_fish_points /
%   load_fish_points_named — only the .X/.Y[/.Z] fields are added.
%
%   PASSTHROUGH: if fish_points(i).pre_transformed is true (e.g. CURVES
%   Format E), the .X and .Y fields are already populated and no rotation
%   is performed.  The struct is returned unchanged.
%
%   METHOD  (follows Transformer.m by Castro-Santos & Goerig 2017, with
%   an added auto-orientation step — see CHANGE NOTE below)
%     For each frame:
%       1. Use the MIDDLE points (all except first and last) to fit a
%          least-squares line in the XY plane: y = a + b*x
%       2. Compute rotation angle:
%            alpha = atan(b)
%            theta = 2*pi - alpha   (clockwise correction)
%       3. Rotate ALL points about the y-intercept (a):
%            x' = x*cos(θ) - (y-a)*sin(θ)
%            y' = x*sin(θ) + (y-a)*cos(θ) + a
%            z' = z  (unchanged — dorso-ventral axis is not rotated)
%       4. AUTO-ORIENTATION (see change note): compare the two ENDPOINTS
%          (point 1 and point end) in rotated space. Whichever has the
%          smaller rotated-X is treated as the head for THIS frame — the
%          code no longer assumes point 1 is always the head. Body length
%          is the distance between the two endpoints, so it is always
%          positive (as long as both endpoints are digitized), instead of
%          silently coming out as exactly zero when a dataset's point
%          order/rotation happens to put the last point ahead of the
%          first.
%            head_ref = min(x'(1), x'(end))
%            tail_ref = max(x'(1), x'(end))
%            bl       = tail_ref - head_ref
%            X = (x' - head_ref) / bl        -> 0 (head) to 1 (tail)
%            Y = (y' - a)        / bl        -> lateral deviation in BL
%            Z = (z  - z_head)   / bl        -> DV deviation in BL (3-D only)
%
%   CHANGE NOTE (why this differs from the original Transformer.m logic):
%   The original method defined "head" as whichever point had the smallest
%   rotated X across ALL digitized points, and body length as the LAST
%   point's position after that shift. That silently breaks in two ways:
%     (a) if a non-endpoint point (e.g. a protruding fin) ends up with a
%         smaller rotated-X than the true head, the origin is wrong; and
%     (b) if the dataset's raw coordinate convention/point order happens
%         to put the LAST point ahead of the FIRST after rotation (as
%         happened with a shark trial where Caudal consistently rotated
%         to a smaller X than Snout), body length comes out as exactly
%         zero on every single frame, and the whole animal's .X/.Y end up
%         100% NaN with no warning.
%   Comparing ONLY the two endpoints (not all points) and taking whichever
%   is smaller as the head fixes both problems and works automatically
%   regardless of point order or coordinate handedness — no need to
%   manually reverse point order in your data.
%
%   Added fields (matrices, rows = frames, cols = points):
%     .X    [nFrames x nPoints]   body-axis position in BL  (0=head, ~1=tail)
%     .Y    [nFrames x nPoints]   lateral displacement in BL (0 = body axis)
%     .Z    [nFrames x nPoints]   dorso-ventral deviation in BL from head (3-D only)
%     .pct_frames_valid   scalar  % of TRACKED frames (rows with any data,
%                                 padding rows excluded) that produced a
%                                 valid transform
%     .n_data_frames      scalar  # of tracked rows (excludes leading/
%                                 trailing all-empty padding rows)
%     .n_frames_reversed  scalar  # of frames where point-end (not point 1)
%                                 was auto-detected as the head — a high
%                                 count here means your point order is
%                                 consistently reversed relative to the
%                                 anatomical head-to-tail convention, which
%                                 is fine (auto-corrected) but worth knowing.

    for fi = 1:numel(fish_points)

        % ---- PASSTHROUGH for pre-transformed data (e.g. CURVES Format E) ----
        if isfield(fish_points(fi), 'pre_transformed') && fish_points(fi).pre_transformed
            fprintf('transform_fish: skipping %s (pre-transformed CURVES data)\n', fish_points(fi).name);
            continue;
        end

        pts     = fish_points(fi).points;   % [nFrames x nPoints x nDims]
        nFrames = size(pts, 1);
        nPoints = size(pts, 2);
        nDims   = size(pts, 3);
        has_z   = (nDims == 3) && isfield(fish_points, 'has_z') && fish_points(fi).has_z;

        % Tracked data window: some CSVs pad the trial with long runs of
        % all-empty rows before/after the actual data. Valid-frame
        % percentages should be computed over rows that contain ANY
        % tracked point, not the padded total row count.
        n_data_frames = sum(any(any(~isnan(pts), 3), 2));

        if nPoints < 3
            error('transform_fish: need at least 3 points; %s has %d.', ...
                  fish_points(fi).name, nPoints);
        end

        middle_idx = 2:nPoints-1;   % exclude head (1) and tail (end)
        % CHANGE NOTE (bug fix): with only 3 tracked points there is a
        % SINGLE middle point and polyfit(x,y,1) silently returns the
        % meaningless min-norm "slope" y/x — the body axis was rotated by
        % an arbitrary angle and body angle was garbage. Fall back to the
        % head-to-tail chord (endpoints) whenever fewer than 2 middle
        % points exist.
        if numel(middle_idx) >= 2
            fit_idx = middle_idx;
        else
            fit_idx = [1, nPoints];
        end

        X = NaN(nFrames, nPoints);
        Y = NaN(nFrames, nPoints);
        Z = NaN(nFrames, nPoints);
        bl_per_frame = NaN(nFrames, 1);
        % CHANGE NOTE (bug fix): "transform_params(nFrames,1) = struct(...)"
        % looks like it preallocates every element with the given NaN
        % defaults, but it does NOT — MATLAB's last-element-assignment
        % growth idiom only fills THAT element; every other element gets
        % [] (empty) in every field, not NaN. That went unnoticed for a
        % while because the main pipeline never reads transform_params
        % (only .X/.Y, which are separately NaN-preallocated and correct).
        % apply_body_transform.m DOES read transform_params directly, and
        % isnan([]) returns [] (not true), so its "skip invalid frames"
        % guard silently failed open on unset frames, leading to an
        % assignment of an empty array into a scalar slot — "Unable to
        % perform assignment because the left and right sides have a
        % different number of elements." repmat() actually copies the
        % NaN-filled struct into every element, unlike last-index growth.
        transform_params = repmat(struct('theta',NaN,'a',NaN,'x1',NaN,'bl',NaN,'sign_flip',NaN), nFrames, 1);

        n_valid    = 0;
        n_reversed = 0;

        for f = 1:nFrames
            x_all = squeeze(pts(f, :, 1));   % [1 x nPoints]
            y_all = squeeze(pts(f, :, 2));
            z_all = [];
            if has_z
                z_all = squeeze(pts(f, :, 3));
            end

            x_mid = x_all(fit_idx);
            y_mid = y_all(fit_idx);

            % Skip frame if any axis-fit point is missing (can't fit the axis line)
            if any(isnan(x_mid)) || any(isnan(y_mid)), continue; end

            % Skip frame if either ENDPOINT is missing — body length is
            % undefined without both ends, regardless of how many middle
            % points are present.
            if isnan(x_all(1)) || isnan(y_all(1)) || isnan(x_all(end)) || isnan(y_all(end))
                continue;
            end

            % Fit line through middle points in XY plane
            coeffs = polyfit(x_mid, y_mid, 1);
            b = coeffs(1);   % slope
            a = coeffs(2);   % y-intercept

            % Rotation angle
            alpha = atan(b);
            theta = 2*pi - alpha;

            % Rotate in XY (Z untouched)
            x_r = x_all .* cos(theta) - (y_all - a) .* sin(theta);
            y_r = x_all .* sin(theta) + (y_all - a) .* cos(theta) + a;

            % ---- Auto-orientation via mirroring (point 1 always = head) ----
            % Body length is the distance between the two endpoints — always
            % non-negative, so this never silently collapses to zero the way
            % "bl = x_trans(end)" could when point-end ended up rotated
            % ahead of point 1 (see CHANGE NOTE above).
            x1 = x_r(1);
            xN = x_r(end);
            bl = abs(xN - x1);

            if isnan(bl) || bl <= eps
                continue;   % degenerate frame (endpoints coincide) — leave NaN
            end

            % Optional known-body-length override (same units as CSV coords)
            if nargin >= 2 && ~isempty(bl_override)
                if isscalar(bl_override)
                    bl = bl_override;
                else
                    bl = bl_override(fi);
                end
                if ~isfinite(bl) || bl <= eps
                    warning(['transform_fish: invalid bl_override for %s — ' ...
                             'ignoring override and using measured length.'], ...
                             fish_points(fi).name);
                    bl = abs(xN - x1);
                end
            end

            % If point-end rotated to a SMALLER x than point 1, mirror the
            % whole frame's X-axis about point 1 so that point 1 still maps
            % to X=0 and point-end still maps to X=1. A pure X-mirror (Y
            % untouched) preserves all pairwise Euclidean distances, so
            % amplitude/curvature/wavelength are unaffected by this — it
            % only fixes which end is labeled "head", not any measured
            % quantity. This keeps point_names{1} correctly meaning "head"
            % for every frame WITHOUT touching your source data's point
            % order or column layout.
            sign_flip = sign(xN - x1);
            if sign_flip == 0
                continue;   % shouldn't happen given bl>eps check above, but guard anyway
            end
            if sign_flip < 0
                n_reversed = n_reversed + 1;
            end

            n_valid = n_valid + 1;

            X(f, :) = sign_flip * (x_r - x1) / bl;   % point1 -> 0, point-end -> 1 always
            Y(f, :) = (y_r - a) / bl;                 % lateral deviation in BL (unaffected by mirror)
            if has_z
                Z(f, :) = (z_all - z_all(1)) / bl;    % DV deviation from point 1 (head) in BL
            end

            bl_per_frame(f) = bl;   % raw body length (original units) — used elsewhere for
                                     % converting real-world displacement to BL/s speed
            transform_params(f).theta     = theta;
            transform_params(f).a         = a;
            transform_params(f).x1        = x1;
            transform_params(f).bl        = bl;
            transform_params(f).sign_flip = sign_flip;
        end

        pct_valid = 100 * n_valid / max(n_data_frames, 1);

        fish_points(fi).X = X;
        fish_points(fi).Y = Y;
        if has_z
            fish_points(fi).Z = Z;
        end
        fish_points(fi).pct_frames_valid  = pct_valid;
        fish_points(fi).n_data_frames     = n_data_frames;     % NEW — tracked rows (excludes padding)
        fish_points(fi).n_frames_reversed = n_reversed;
        fish_points(fi).bl_per_frame      = bl_per_frame;      % NEW — raw body length per frame
        fish_points(fi).transform_params  = transform_params;  % NEW — for projecting other points (fin roots, girdle pts) into the same body-relative frame

        fprintf('transform_fish: %s | %d/%d tracked frames valid (%.1f%%)', ...
                fish_points(fi).name, n_valid, n_data_frames, pct_valid);
        if n_reversed > 0
            fprintf('  [%d/%d frame(s) mirrored: point-end rotated ahead of point 1]', ...
                    n_reversed, max(n_valid,1));
        end
        fprintf('\n');
        if n_valid > 0 && n_reversed > n_valid / 2
            fprintf(['  NOTE: point-end was rotated ahead of point 1 in the majority of frames\n' ...
                     '  for %s. Each such frame was auto-mirrored so point_names{1} ("%s") still\n' ...
                     '  correctly maps to the head (X~0) and point_names{end} ("%s") to the tail\n' ...
                     '  (X~1) — amplitude/curvature are unaffected by this (distances are preserved\n' ...
                     '  under mirroring). No changes to your source CSV were needed.\n'], ...
                     fish_points(fi).name, fish_points(fi).point_names{1}, fish_points(fi).point_names{end});
        end

        if n_valid == 0
            warning(['transform_fish: %s has ZERO valid frames after transform. ' ...
                     'All downstream kinematics for this animal will be NaN — ' ...
                     'check that middle points and both endpoints are actually ' ...
                     'tracked in at least some frames.'], fish_points(fi).name);
        elseif pct_valid < 20
            warning(['transform_fish: %s has only %.1f%% valid frames — results ' ...
                     'may be based on very little data. Check tracking coverage.'], ...
                     fish_points(fi).name, pct_valid);
        end
    end
end
