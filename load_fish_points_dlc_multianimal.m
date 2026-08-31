function fish_points = load_fish_points_dlc_multianimal(filename, p_cutoff)
% LOAD_FISH_POINTS_DLC_MULTIANIMAL  Load a native DeepLabCut multi-animal
%                                    project CSV export.
%
%   fish_points = load_fish_points_dlc_multianimal(filename)
%   fish_points = load_fish_points_dlc_multianimal(filename, p_cutoff)
%
%   This is DIFFERENT from every other loader in this toolkit: a DLC
%   multi-animal project export has FOUR header rows, not one or three:
%       row1: scorer,       DLC_resnet50_..., DLC_resnet50_..., ...
%       row2: individuals,  individual1,      individual1,     ...
%       row3: bodyparts,    snout,            snout,           ...
%       row4: coords,       x,                y,          likelihood, ...
%       row5+: <frame idx>, <data...>
%   Every fish's worth of columns repeats for each unique value in the
%   "individuals" row. Do NOT try load_fish_points_named on this format;
%   its single-row-header assumption doesn't match a 4-row DLC header at
%   all (same class of mismatch as the raw single-animal 3-row DLC header
%   handled elsewhere in this toolkit; this format just has one MORE
%   header row on top of that).
%
%   INPUTS
%     filename   - path to the DLC multi-animal CSV
%     p_cutoff   - (optional) likelihood threshold in [0,1]; any x/y with
%                  likelihood below this is set to NaN. Default 0 (no
%                  filtering).
%
%   OUTPUT  fish_points: struct ARRAY, one element per unique individual
%   found in the "individuals" row, in the SAME schema as
%   load_fish_points.m, a drop-in for anything expecting that struct
%   array, including compute_polarization / compute_angle_to_flow /
%   compute_distance_between_individuals / filter_dlc_jumps:
%     .name:       the individual's label (e.g. 'individual1')
%     .frames      [nFrames x 1]
%     .point_names {1 x nPoints}  bodypart names, in first-appearance
%                  order for that individual (e.g. 'snout','mid 1',...)
%     .points      [nFrames x nPoints x 2]  (x,y), this format doesn't
%                  carry a Z coordinate
%     .has_z       false
%     .format      'dlc_multianimal'

    if nargin < 2 || isempty(p_cutoff), p_cutoff = 0; end

    raw = readcell(filename);   % preserves the mixed text/number header rows
    row_individuals = raw(2, :);
    row_bodyparts   = raw(3, :);
    row_coords      = raw(4, :);
    data            = raw(5:end, :);   % column numbering UNCHANGED from raw (col 1 = frame idx)
    nRows = size(data, 1);

    frames = cellfun(@to_num, data(:,1));

    % Column 1 is the frame index in every row above, individuals/
    % bodyparts/coords are only meaningful from column 2 onward.
    dataCols = 2:size(raw,2);
    individuals = row_individuals(dataCols);
    bodyparts   = row_bodyparts(dataCols);
    coords      = row_coords(dataCols);

    uniq_individuals = unique(individuals, 'stable');
    nFish = numel(uniq_individuals);

    if nFish == 0
        error('load_fish_points_dlc_multianimal: no individuals found in row 2 of %s: check this is really a multi-animal DLC export.', filename);
    end

    fish_points(nFish) = struct('name','', 'frames',[], 'point_names',{{}}, ...
                                  'points',[], 'has_z',false, 'format','');

    for fi = 1:nFish
        indiv = uniq_individuals{fi};
        cols_this_fish = find(strcmp(individuals, indiv));   % indices INTO dataCols (i.e. relative to column 2)

        bp_this = bodyparts(cols_this_fish);
        point_names = unique(bp_this, 'stable');
        nPoints = numel(point_names);

        pts = NaN(nRows, nPoints, 2);
        n_missing = 0;

        for pi = 1:nPoints
            bp = point_names{pi};
            match = cols_this_fish(strcmp(bodyparts(cols_this_fish), bp));

            col_x = match(strcmpi(coords(match), 'x'));
            col_y = match(strcmpi(coords(match), 'y'));
            col_p = match(strcmpi(coords(match), 'likelihood'));

            if isempty(col_x) || isempty(col_y)
                warning('load_fish_points_dlc_multianimal: bodypart "%s" missing x/y for %s: left as NaN.', bp, indiv);
                n_missing = n_missing + 1;
                continue;
            end

            % +1 to convert from "relative to column 2" back to actual
            % column number in `data` (which retains raw's numbering,
            % including column 1 = frame index).
            x = cellfun(@to_num, data(:, col_x(1)+1));
            y = cellfun(@to_num, data(:, col_y(1)+1));

            if ~isempty(col_p) && p_cutoff > 0
                p = cellfun(@to_num, data(:, col_p(1)+1));
                low_conf = p < p_cutoff;
                x(low_conf) = NaN;
                y(low_conf) = NaN;
            end

            pts(:, pi, 1) = x;
            pts(:, pi, 2) = y;
        end

        fish_points(fi).name        = indiv;
        fish_points(fi).frames      = frames;
        fish_points(fi).point_names = point_names;
        fish_points(fi).points      = pts;
        fish_points(fi).has_z       = false;
        fish_points(fi).format      = 'dlc_multianimal';

        n_valid = sum(~any(isnan(pts(:,:,1)) | isnan(pts(:,:,2)), 2));
        fprintf('Loaded (DLC multi-animal): %s  [2D, %d frames, %d points, p-cutoff=%.2f, %d/%d frames fully valid]\n', ...
                indiv, nRows, nPoints, p_cutoff, n_valid, nRows);
        fprintf('  Points: %s\n', strjoin(point_names, ' | '));
    end

    fprintf('load_fish_points_dlc_multianimal: %s -> %d individual(s) found\n', filename, nFish);
end


function v = to_num(x)
    if isnumeric(x)
        v = x;
    elseif ischar(x) || isstring(x)
        v = str2double(x);
    else
        v = NaN;
    end
end
