function fish_points = load_fish_points_named(filename, selected_points, point_order)
% LOAD_FISH_POINTS_NAMED  Load a named-column tracking CSV (e.g. eye_X, snout_Y, snout_Z).
%
%   fish_points = load_fish_points_named(filename, selected_points, point_order)
%   fish_points = load_fish_points_named(filename)   % returns available point names only
%
%   INPUTS
%     filename         - path to CSV.  Columns follow the pattern <name>_X, <name>_Y[, <name>_Z].
%                        One file = one animal (no fish-grouping column).
%                        Row index is used as frame number.
%
%     selected_points  - cell array of point base-names to include, e.g.:
%                          {'snout', 'peduncle', 'caudaltip'}
%                        Each entry may also be a two-element cell meaning
%                        "average of these two":
%                          {'snout', {'Rpectbase','Lpectbase'}, 'peduncle', 'caudaltip'}
%                        If omitted or empty, the function prints available names and returns.
%
%     point_order      - integer vector giving the order of selected_points
%                        along the body (1 = head, end = tail).  E.g. [1 2 3 4].
%                        If omitted, the order of selected_points is used as-is.
%
%   OUTPUT  fish_points — 1-element struct (one animal per file):
%     .name        string  — filename stem
%     .frames      [nFrames x 1]
%     .point_names {1 x nPoints}  ordered labels (averaged pairs shown as 'A+B')
%     .points      [nFrames x nPoints x nDims]   nDims = 2 or 3
%     .has_z       logical
%
%   EXAMPLE
%     % Step 1: discover available points
%     load_fish_points_named('poly1_data.csv');
%
%     % Step 2: load with selection and ordering
%     fp = load_fish_points_named('poly1_data.csv', ...
%            {'snout', {'Rpectbase','Lpectbase'}, 'peduncle', 'caudaltip'}, ...
%            [1 2 3 4]);

    %% Read file
    opts = detectImportOptions(filename);
    opts.VariableNamingRule = 'preserve';
    T = readtable(filename, opts);
    colNames = T.Properties.VariableNames;
    nFrames  = height(T);
    frames   = (1:nFrames)';

    %% Detect available point names and dimensionality
    % Columns: <name>_X, <name>_Y, <name>_Z  (case-insensitive suffix)
    tok_x = regexp(colNames, '^(.+)_[Xx]$', 'tokens');
    tok_z = regexp(colNames, '^(.+)_[Zz]$', 'tokens');
    has_z = any(~cellfun(@isempty, tok_z));

    base_names = cellfun(@(t) t{1}{1}, tok_x(~cellfun(@isempty,tok_x)), 'UniformOutput', false);
    base_names = unique(base_names, 'stable');

    %% If no selection given — print available names and return
    if nargin < 2 || isempty(selected_points)
        fprintf('\nAvailable points in %s:\n', filename);
        for i = 1:numel(base_names)
            fprintf('  %2d.  %s\n', i, base_names{i});
        end
        fprintf('\nHas Z dimension: %s\n\n', mat2str(has_z));
        fish_points = struct('name', filename, 'frames', frames, ...
                             'point_names', {base_names}, 'points', [], 'has_z', has_z);
        return
    end

    %% Apply ordering
    if nargin < 3 || isempty(point_order)
        point_order = 1:numel(selected_points);
    end
    selected_points = selected_points(point_order);

    %% Build points array
    nDims   = 2 + has_z;
    nPoints = numel(selected_points);
    pts     = NaN(nFrames, nPoints, nDims);
    labels  = cell(1, nPoints);
    dims    = {'X','Y','Z'};

    for pi = 1:nPoints
        entry = selected_points{pi};

        if iscell(entry)
            % Average of two (or more) landmarks
            pair_data = NaN(nFrames, numel(entry), nDims);
            for ei = 1:numel(entry)
                for di = 1:nDims
                    col = find_col(colNames, entry{ei}, dims{di});
                    if ~isempty(col)
                        pair_data(:, ei, di) = T.(colNames{col});
                    end
                end
            end
            pts(:, pi, :) = mean(pair_data, 2, 'omitnan');
            labels{pi}    = strjoin(entry, '+');
        else
            % Single landmark
            for di = 1:nDims
                col = find_col(colNames, entry, dims{di});
                if ~isempty(col)
                    pts(:, pi, di) = T.(colNames{col});
                end
            end
            labels{pi} = entry;
        end
    end

    %% Build output struct (same schema as load_fish_points)
    [~, stem] = fileparts(filename);
    fish_points = struct( ...
        'name',        stem, ...
        'frames',      frames, ...
        'point_names', {labels}, ...
        'points',      pts, ...
        'has_z',       has_z);

    nDimStr = sprintf('%dD', nDims);
    fprintf('Loaded: %s  [%s, %d frames, %d points]\n', stem, nDimStr, nFrames, nPoints);
    fprintf('  Points: %s\n', strjoin(labels, ' → '));
end

% -------------------------------------------------------------------------
function idx = find_col(colNames, base, dim)
% Case-insensitive search for <base>_<dim> in colNames.
    pattern = sprintf('^%s_%s$', regexptranslate('escape', base), dim);
    idx = find(~cellfun(@isempty, regexpi(colNames, pattern)), 1);
end
