function dist = compute_distance_between_individuals(fish_points, snoutName, cm_per_unit)
% COMPUTE_DISTANCE_BETWEEN_INDIVIDUALS  Pairwise snout-to-snout distance
%                                       between all fish, per frame.
%
%   dist = compute_distance_between_individuals(fish_points, snoutName)
%   dist = compute_distance_between_individuals(fish_points, snoutName, cm_per_unit)
%
%   Pass the RAW (untransformed) fish_points struct array — same
%   reasoning as compute_polarization.m: this needs everyone in the same
%   shared (camera/world) coordinate frame, not each fish's own private
%   post-transform_fish frame.
%
%   INPUTS
%     fish_points   - struct ARRAY, one element per fish, RAW
%                     (pre-transform_fish). Each needs .points [nFrames x
%                     nPoints x nDims>=2] and .point_names.
%     snoutName     - point name for the snout landmark
%     cm_per_unit   - (optional) conversion factor from your CSV's raw
%                     coordinate units to cm (e.g. if your DLC output is
%                     already calibrated to cm, leave this at the default
%                     of 1.0; if it's in pixels, this needs to be your
%                     pixel-to-cm calibration factor for the distances
%                     below to actually be in cm as intended). Default 1.0
%                     — i.e. distances are reported in whatever units your
%                     raw coordinates are already in unless you supply this.
%
%   OUTPUT  dist — struct:
%     .pairwise_dist       [nFrames x nFish x nFish] symmetric distance
%                          matrix per frame (NaN on the diagonal and
%                          wherever a fish's snout isn't tracked that frame)
%     .fish_names          {1 x nFish} cell array of names, matching the
%                          pairwise_dist indexing
%     .mean_nn_dist        [nFrames x 1] mean NEAREST-NEIGHBOR distance
%                          per frame — each fish's distance to its single
%                          closest neighbor, averaged across fish present
%     .mean_pairwise_dist  [nFrames x 1] mean of ALL pairwise distances
%                          per frame (not just nearest-neighbor)
%     .overall_mean_nn_dist / .overall_mean_pairwise_dist
%                          summary over all valid frames
%     .cm_per_unit         the conversion factor actually used

    if nargin < 3 || isempty(cm_per_unit), cm_per_unit = 1.0; end

    nFish = numel(fish_points);
    if nFish < 2
        warning(['compute_distance_between_individuals: only %d fish provided — ' ...
                 'inter-individual distance needs at least 2.'], nFish);
    end

    nFrames = size(fish_points(1).points, 1);
    snout_xy = NaN(nFrames, nFish, 2);
    fish_names = cell(1, nFish);

    for k = 1:nFish
        fish_names{k} = safe_name(fish_points(k));
        pn = fish_points(k).point_names;
        iSnout = find(strcmpi(pn, snoutName), 1);
        if isempty(iSnout)
            warning('compute_distance_between_individuals: fish %d ("%s") missing "%s" — excluded.', ...
                     k, fish_names{k}, snoutName);
            continue;
        end
        nF_k = min(nFrames, size(fish_points(k).points,1));
        snout_xy(1:nF_k, k, :) = fish_points(k).points(1:nF_k, iSnout, 1:2);
    end

    pairwise_dist = NaN(nFrames, nFish, nFish);
    mean_nn_dist = NaN(nFrames, 1);
    mean_pairwise_dist = NaN(nFrames, 1);

    for f = 1:nFrames
        pos = squeeze(snout_xy(f, :, :));   % [nFish x 2]
        if nFish == 1, pos = pos(:)'; end   % squeeze quirk guard for nFish==1
        D = NaN(nFish, nFish);
        for i = 1:nFish
            for j = 1:nFish
                if i == j, continue; end
                if any(isnan(pos(i,:))) || any(isnan(pos(j,:))), continue; end
                D(i,j) = norm(pos(i,:) - pos(j,:)) * cm_per_unit;
            end
        end
        pairwise_dist(f,:,:) = D;

        nn = min(D, [], 2, 'omitnan');   % each fish's nearest-neighbor distance
        nn = nn(~isinf(nn) & ~isnan(nn));
        if ~isempty(nn), mean_nn_dist(f) = mean(nn); end

        upper = D(triu(true(nFish),1));
        upper = upper(~isnan(upper));
        if ~isempty(upper), mean_pairwise_dist(f) = mean(upper); end
    end

    dist.pairwise_dist      = pairwise_dist;
    dist.fish_names         = fish_names;
    dist.mean_nn_dist        = mean_nn_dist;
    dist.mean_pairwise_dist  = mean_pairwise_dist;
    dist.overall_mean_nn_dist       = mean(mean_nn_dist, 'omitnan');
    dist.overall_mean_pairwise_dist = mean(mean_pairwise_dist, 'omitnan');
    dist.cm_per_unit         = cm_per_unit;

    fprintf(['Distance between individuals: %d fish | mean nearest-neighbor dist = %s, ' ...
             'mean pairwise dist = %s (units: raw * %.4g)\n'], nFish, ...
            fmt_val(dist.overall_mean_nn_dist), fmt_val(dist.overall_mean_pairwise_dist), cm_per_unit);
    if cm_per_unit == 1.0
        fprintf(['  NOTE: cm_per_unit was left at the default (1.0) — these distances are in ' ...
                 'whatever units your raw coordinates are in, NOT necessarily cm. Pass your ' ...
                 'pixel-to-cm (or other unit-to-cm) calibration factor if you need real cm.\n']);
    end
end

function s = safe_name(fp)
    if isfield(fp, 'name') && ~isempty(fp.name), s = fp.name; else, s = 'unnamed'; end
end

function s = fmt_val(v)
    if isnan(v), s = 'NaN'; else, s = sprintf('%.4f', v); end
end
