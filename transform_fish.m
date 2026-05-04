function fish_points = transform_fish(fish_points)
% TRANSFORM_FISH  Rotate and translate each fish's midline so the mean
%                 swimming axis is parallel to the X-axis, with the first
%                 point (head) at X = 0.
%
%   fish_points = transform_fish(fish_points)
%
%   Works with any number of points (>= 3) and with 2-D or 3-D data.
%   The struct schema is identical to the output of load_fish_points /
%   load_fish_points_named — only the .X/.Y[/.Z] fields are added.
%
%   METHOD  (follows Transformer.m by Castro-Santos & Goerig 2017)
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
%       4. Translate so the minimum x' (head) = 0:
%            X = x' - min(x')
%            Y = y'
%            Z = z'  (if present)
%
%   Added fields (matrices, rows = frames, cols = points):
%     .X    [nFrames x nPoints]   rotated + translated x
%     .Y    [nFrames x nPoints]   rotated y
%     .Z    [nFrames x nPoints]   z (3-D only, unchanged by rotation)

    for fi = 1:numel(fish_points)
        pts     = fish_points(fi).points;   % [nFrames x nPoints x nDims]
        nFrames = size(pts, 1);
        nPoints = size(pts, 2);
        nDims   = size(pts, 3);
        has_z   = (nDims == 3) && isfield(fish_points, 'has_z') && fish_points(fi).has_z;

        if nPoints < 3
            error('transform_fish: need at least 3 points; %s has %d.', ...
                  fish_points(fi).name, nPoints);
        end

        middle_idx = 2:nPoints-1;   % exclude head (1) and tail (end)

        X = NaN(nFrames, nPoints);
        Y = NaN(nFrames, nPoints);
        Z = NaN(nFrames, nPoints);

        for f = 1:nFrames
            x_all = squeeze(pts(f, :, 1));   % [1 x nPoints]
            y_all = squeeze(pts(f, :, 2));
            z_all = [];
            if has_z
                z_all = squeeze(pts(f, :, 3));
            end

            x_mid = x_all(middle_idx);
            y_mid = y_all(middle_idx);

            % Skip frame if any middle point is missing
            if any(isnan(x_mid)) || any(isnan(y_mid)), continue; end

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

            % Translate: head (min x) to x = 0
            x_shift = min(x_r);
            X(f, :) = x_r - x_shift;
            Y(f, :) = y_r;
            if has_z
                Z(f, :) = z_all;   % Z not rotated
            end
        end

        fish_points(fi).X = X;
        fish_points(fi).Y = Y;
        if has_z
            fish_points(fi).Z = Z;
        end
    end
end
