function pol = compute_polarization(fish_points, snoutName, peduncleName, flowAxisDeg)
% COMPUTE_POLARIZATION  School-level polarization (heading alignment)
%                       from snout+peduncle vectors across multiple fish.
%
%   pol = compute_polarization(fish_points, snoutName, peduncleName)
%   pol = compute_polarization(fish_points, snoutName, peduncleName, flowAxisDeg)
%
%   IMPORTANT: pass the RAW (untransformed) fish_points struct array —
%   i.e. straight from load_fish_points()/load_fish_points_named(),
%   BEFORE transform_fish(). transform_fish deliberately puts each fish
%   into its OWN private body-length-normalized coordinate frame, which
%   destroys the shared spatial/angular relationships between fish that
%   polarization needs. This function needs everyone in the SAME
%   (camera/world) coordinate frame.
%
%   INPUTS
%     fish_points   - struct ARRAY, one element per fish (as returned by
%                     load_fish_points() for multi-animal Format A, or
%                     assembled by hand from multiple single-fish loads —
%                     see note below). Each element needs .points
%                     [nFrames x nPoints x nDims>=2] and .point_names.
%     snoutName     - point name (or numbered label) for the snout landmark
%     peduncleName  - point name (or numbered label) for the peduncle landmark
%     flowAxisDeg   - (optional) the heading angle (degrees) that counts
%                     as "aligned with flow," e.g. 0 for left-to-right
%                     flow in raw image/camera coordinates. Default 0.
%                     This only affects .heading_to_flow_deg (per-fish);
%                     polarization itself doesn't depend on flow
%                     direction — it measures alignment of fish WITH
%                     EACH OTHER, regardless of which way that is.
%
%   OUTPUT  pol — struct:
%     .heading_deg          [nFrames x nFish] each fish's heading angle
%                            (peduncle -> snout vector), degrees, in the
%                            raw/world coordinate frame
%     .heading_to_flow_deg  [nFrames x nFish] heading_deg re-centered so
%                            0 = aligned with flowAxisDeg
%     .polarization         [nFrames x 1] school-level alignment metric,
%                            0 (headings random/opposed) to 1 (all fish
%                            facing exactly the same way). NaN in frames
%                            with fewer than 2 fish present.
%     .n_fish_present       [nFrames x 1] how many fish had both snout
%                            and peduncle tracked that frame
%     .mean_polarization / .std_polarization   summary over valid frames
%
%   METHOD: standard circular-statistics polarization —
%     polarization(f) = | mean_k( exp(i * heading_k(f)) ) |
%   i.e. the length of the average unit heading vector across all fish
%   present in frame f. This is invariant to the OVERALL direction
%   they're all facing — it only measures how well-aligned they are WITH
%   EACH OTHER.
%
%   NOTE ON MULTI-FISH DATA ASSEMBLY: if your fish are in separate CSVs
%   or you loaded them individually, build the struct array by hand
%   before calling this, e.g.:
%       fp(1) = load_fish_points_named('fish1.csv', {'snout','peduncle'}, [1 2]);
%       fp(2) = load_fish_points_named('fish2.csv', {'snout','peduncle'}, [1 2]);
%       pol = compute_polarization(fp, 'snout', 'peduncle');
%   All fish must share the same frame indexing (frame N in fish 1 is the
%   same real time as frame N in fish 2) for this to be meaningful.

    if nargin < 4 || isempty(flowAxisDeg), flowAxisDeg = 0; end

    nFish = numel(fish_points);
    if nFish < 2
        warning(['compute_polarization: only %d fish provided — polarization is a ' ...
                 'group-alignment metric and needs at least 2 to mean anything.'], nFish);
    end

    nFrames = size(fish_points(1).points, 1);
    heading_deg = NaN(nFrames, nFish);

    for k = 1:nFish
        pn = fish_points(k).point_names;
        iSnout = find(strcmpi(pn, snoutName), 1);
        iPed   = find(strcmpi(pn, peduncleName), 1);
        if isempty(iSnout) || isempty(iPed)
            warning(['compute_polarization: fish %d ("%s") is missing "%s" and/or "%s" — ' ...
                     'excluded from polarization for all frames.'], ...
                     k, fish_points(k).name, snoutName, peduncleName);
            continue;
        end
        if size(fish_points(k).points,1) ~= nFrames
            warning(['compute_polarization: fish %d ("%s") has %d frames, expected %d ' ...
                     '(from fish 1) — frame indices may not correspond to the same real ' ...
                     'time across fish. Proceeding, but verify your inputs are synchronized.'], ...
                     k, fish_points(k).name, size(fish_points(k).points,1), nFrames);
        end

        nF_k = min(nFrames, size(fish_points(k).points,1));
        snout_xy = squeeze(fish_points(k).points(1:nF_k, iSnout, 1:2));
        ped_xy   = squeeze(fish_points(k).points(1:nF_k, iPed,   1:2));

        vx = snout_xy(:,1) - ped_xy(:,1);
        vy = snout_xy(:,2) - ped_xy(:,2);
        heading_deg(1:nF_k, k) = atan2d(vy, vx);
    end

    % ---- Per-frame polarization ----
    n_fish_present = sum(~isnan(heading_deg), 2);
    polarization = NaN(nFrames, 1);
    for f = 1:nFrames
        theta = heading_deg(f, ~isnan(heading_deg(f,:)));
        if numel(theta) >= 2
            polarization(f) = abs(mean(cosd(theta)) + 1i*mean(sind(theta)));
        end
    end

    valid = ~isnan(polarization);
    if any(valid)
        mean_pol = mean(polarization(valid));
        std_pol  = std(polarization(valid));
    else
        mean_pol = NaN; std_pol = NaN;
    end

    heading_to_flow = mod(heading_deg - flowAxisDeg + 180, 360) - 180;   % wrap to (-180,180]

    pol.heading_deg         = heading_deg;
    pol.heading_to_flow_deg = heading_to_flow;
    pol.polarization        = polarization;
    pol.n_fish_present      = n_fish_present;
    pol.mean_polarization   = mean_pol;
    pol.std_polarization    = std_pol;
    pol.flowAxisDeg         = flowAxisDeg;
    pol.n_fish              = nFish;

    fprintf('Polarization: %d fish | mean polarization = %s (SD %s) over %d/%d frames with >=2 fish present\n', ...
            nFish, fmt_val(mean_pol), fmt_val(std_pol), sum(n_fish_present>=2), nFrames);
end

function s = fmt_val(v)
    if isnan(v), s = 'NaN'; else, s = sprintf('%.4f', v); end
end
