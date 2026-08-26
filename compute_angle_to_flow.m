function ang = compute_angle_to_flow(fish_points, snoutName, peduncleName, flowAxisDeg)
% COMPUTE_ANGLE_TO_FLOW  Per-fish heading angle relative to a flow direction.
%
%   ang = compute_angle_to_flow(fish_points, snoutName, peduncleName)
%   ang = compute_angle_to_flow(fish_points, snoutName, peduncleName, flowAxisDeg)
%
%   Works for one fish or a struct array of several — each is handled
%   independently (this is NOT a group metric; see compute_polarization.m
%   for that). Pass the RAW (untransformed) fish_points — same reasoning
%   as compute_polarization.m: transform_fish's per-fish normalized frame
%   would distort the angle relative to a shared flow direction.
%
%   INPUTS
%     fish_points   - struct or struct array, RAW (pre-transform_fish).
%                     Each element needs .points [nFrames x nPoints x
%                     nDims>=2] and .point_names.
%     snoutName     - point name for the snout landmark
%     peduncleName  - point name for the peduncle landmark
%     flowAxisDeg   - (optional) heading angle (deg) that counts as
%                     "aligned with flow" in the raw/camera coordinate
%                     frame, e.g. 0 for left-to-right flow. Default 0.
%
%   OUTPUT  ang — struct ARRAY (one per fish):
%     .name                  fish name (from fish_points(k).name if present)
%     .heading_deg           [nFrames x 1] raw heading angle (peduncle->snout)
%     .angle_to_flow_deg     [nFrames x 1] heading re-centered so 0 = aligned
%                            with flow, wrapped to (-180, 180]. Positive/
%                            negative sign follows standard atan2d
%                            convention (counter-clockwise positive) —
%                            check against your own footage's orientation
%                            if up/down or left/right sign matters for
%                            your analysis.
%     .mean_angle_to_flow_deg / .std_angle_to_flow_deg / .range_angle_to_flow_deg
%         CIRCULAR statistics — the mean is the direction of the average
%         unit heading vector (atan2d of the summed sines/cosines) and the
%         std is the circular standard deviation sqrt(-2*ln R) in deg,
%         where R is the mean resultant length. A plain arithmetic mean/std
%         of wrapped angles is garbage whenever headings straddle the
%         +/-180 wrap (e.g. fish swimming straight into the flow at ~180
%         deg: +179/-179 frame noise averages to ~0 with an SD of ~180,
%         i.e. "pointing in every direction" while polarization says
%         otherwise). Range is the total angular span of the data measured
%         relative to the circular mean, so it is also wrap-safe.
%     .n_valid_frames

    if nargin < 4 || isempty(flowAxisDeg), flowAxisDeg = 0; end

    nFish = numel(fish_points);
    ang(nFish) = struct();

    for k = 1:nFish
        pn = fish_points(k).point_names;
        iSnout = find(strcmpi(pn, snoutName), 1);
        iPed   = find(strcmpi(pn, peduncleName), 1);

        nFrames = size(fish_points(k).points, 1);
        heading_deg = NaN(nFrames, 1);

        if isempty(iSnout) || isempty(iPed)
            warning(['compute_angle_to_flow: fish %d ("%s") is missing "%s" and/or "%s" — ' ...
                     'all frames NaN for this fish.'], k, safe_name(fish_points(k)), snoutName, peduncleName);
        else
            snout_xy = squeeze(fish_points(k).points(:, iSnout, 1:2));
            ped_xy   = squeeze(fish_points(k).points(:, iPed,   1:2));
            vx = snout_xy(:,1) - ped_xy(:,1);
            vy = snout_xy(:,2) - ped_xy(:,2);
            heading_deg = atan2d(vy, vx);
        end

        angle_to_flow = mod(heading_deg - flowAxisDeg + 180, 360) - 180;

        valid = ~isnan(angle_to_flow);
        if any(valid)
            a = angle_to_flow(valid);
            % Circular mean = direction of the average unit heading vector
            mean_a = atan2d(sum(sind(a)), sum(cosd(a)));
            R = hypot(sum(cosd(a)), sum(sind(a))) / numel(a);   % resultant length, 0..1
            if R > 0
                std_a = sqrt(-2*log(R)) * 180/pi;   % circular SD in deg
            else
                std_a = 180;   % perfectly scattered headings
            end
            % Wrap-safe range: angular span of the data around the
            % circular mean (so 179/-179 gives ~2 deg, not 358).
            dev = mod(a - mean_a + 180, 360) - 180;
            range_a = max(dev) - min(dev);
        else
            mean_a = NaN; std_a = NaN; range_a = NaN;
        end

        ang(k).name                    = safe_name(fish_points(k));
        ang(k).heading_deg             = heading_deg;
        ang(k).angle_to_flow_deg       = angle_to_flow;
        ang(k).mean_angle_to_flow_deg  = mean_a;
        ang(k).std_angle_to_flow_deg   = std_a;
        ang(k).range_angle_to_flow_deg = range_a;
        ang(k).n_valid_frames          = sum(valid);

        fprintf('Angle to flow (%s): mean=%s\xB1%s deg (range %s), %d/%d valid frames\n', ...
                ang(k).name, fmt_val(mean_a), fmt_val(std_a), fmt_val(range_a), sum(valid), nFrames);
    end
end

function s = safe_name(fp)
    if isfield(fp, 'name') && ~isempty(fp.name), s = fp.name; else, s = 'unnamed'; end
end

function s = fmt_val(v)
    if isnan(v), s = 'NaN'; else, s = sprintf('%.2f', v); end
end
