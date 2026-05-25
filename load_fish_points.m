function fish_points = load_fish_points(filename)
% LOAD_FISH_POINTS  Load a DLC-style or numbered-point tracking CSV.
%
%   fish_points = load_fish_points(filename)
%
%   Supported column formats:
%     FORMAT A (DLC-style):      Fish1_P1_x, Fish1_P1_y[, Fish1_P1_z]
%     FORMAT C (numbered 3D/2D): pt1_X, pt1_Y[, pt1_Z], pt2_X, pt2_Y[, pt2_Z], ...
%
%   For FORMAT C, the predefined point anatomy is:
%     pt1  = base of pectoral fin (right)
%     pt2  = tip of pectoral fin (right)
%     pt3  = peduncle
%     pt4  = tip of caudal fin
%     pt5  = eye (right)
%     pt6  = base of pelvic fin
%     pt7  = tip of pelvic fin
%     pt8  = base of anal fin
%     pt9  = tip of anal fin
%     pt10 = tip of dorsal fin
%     pt11 = eye (left / #2)
%     pt12 = pectoral fin #2 (left base) — 2D only
%
%   OUTPUT  fish_points — struct array, one element per fish:
%     .name        string
%     .frames      [nFrames x 1]
%     .point_names {1 x nPoints}  cell of point labels
%     .points      [nFrames x nPoints x nDims]   nDims = 2 or 3
%     .has_z       logical scalar
%     .format      string: 'DLC' or 'numbered'

    opts = detectImportOptions(filename);
    opts.VariableNamingRule = 'preserve';
    T = readtable(filename, opts);

    colNames = T.Properties.VariableNames;
    nFrames  = height(T);

    %% ---- Detect format ----

    % FORMAT C: pt<N>_X  (numbered points, no fish grouping)
    tok_ptnum = regexp(colNames, '^pt(\d+)_[XxYyZz]$', 'tokens');
    is_ptnum  = ~cellfun(@isempty, tok_ptnum);

    if any(is_ptnum)
        fish_points = load_numbered_pts(T, colNames, nFrames, filename);
        return;
    end

    %% ---- FORMAT A: DLC-style Fish1_P1_x ----
    frames = T.frame;

    tok_xy  = regexp(colNames, '^(Fish\d+)_(P\d+)_([xy])$',  'tokens');
    tok_xyz = regexp(colNames, '^(Fish\d+)_(P\d+)_([xyz])$', 'tokens');
    has_z   = any(~cellfun(@isempty, tok_xyz) & cellfun(@isempty, tok_xy));

    tok = tok_xyz;
    valid = ~cellfun(@isempty, tok);

    fish_names  = unique(cellfun(@(t) t{1}{1}, tok(valid), 'UniformOutput', false), 'stable');
    point_names = unique(cellfun(@(t) t{1}{2}, tok(valid), 'UniformOutput', false), 'stable');
    nFish   = numel(fish_names);
    nPoints = numel(point_names);
    nDims   = 2 + has_z;

    fish_points(nFish) = struct('name','','frames',[],'point_names',{{}},'points',[],'has_z',false,'format','');

    for fi = 1:nFish
        fname = fish_names{fi};
        pts   = NaN(nFrames, nPoints, nDims);
        dims  = {'x','y','z'};

        for pi = 1:nPoints
            for di = 1:nDims
                col = sprintf('%s_%s_%s', fname, point_names{pi}, dims{di});
                if ismember(col, colNames)
                    pts(:, pi, di) = T.(col);
                end
            end
        end

        fish_points(fi).name        = fname;
        fish_points(fi).frames      = frames;
        fish_points(fi).point_names = point_names;
        fish_points(fi).points      = pts;
        fish_points(fi).has_z       = has_z;
        fish_points(fi).format      = 'DLC';
    end

    nDimStr = sprintf('%dD', nDims);
    fprintf('Loaded: %s  [%s, %d frames, %d fish, %d points]  [DLC format]\n', ...
            filename, nDimStr, nFrames, nFish, nPoints);
end


% =========================================================================
%  FORMAT C LOADER: pt1_X, pt1_Y[, pt1_Z], pt2_X, ...
% =========================================================================
function fish_points = load_numbered_pts(T, colNames, nFrames, filename)
% Load numbered-point format and attach anatomical names.

    POINT_ANATOMY = { ...
        'pt1',  'pect_base_R'; ...
        'pt2',  'pect_tip_R'; ...
        'pt3',  'peduncle'; ...
        'pt4',  'caudal_tip'; ...
        'pt5',  'eye_R'; ...
        'pt6',  'pelvic_base'; ...
        'pt7',  'pelvic_tip'; ...
        'pt8',  'anal_base'; ...
        'pt9',  'anal_tip'; ...
        'pt10', 'dorsal_tip'; ...
        'pt11', 'eye_L'; ...
        'pt12', 'pect_base_L'; ...
    };

    % Find all unique point numbers present
    tok_all = regexp(colNames, '^(pt\d+)_[XxYyZz]$', 'tokens');
    pt_bases = cellfun(@(t) t{1}{1}, tok_all(~cellfun(@isempty, tok_all)), 'UniformOutput', false);
    pt_bases = unique(pt_bases, 'stable');

    % Sort by numeric order
    pt_nums  = cellfun(@(s) str2double(regexp(s,'\d+','match','once')), pt_bases);
    [~, ord] = sort(pt_nums);
    pt_bases = pt_bases(ord);
    nPoints  = numel(pt_bases);

    % Determine if Z exists
    has_z = any(~cellfun(@isempty, regexp(colNames, '^pt\d+_[Zz]$')));
    nDims = 2 + has_z;

    pts    = NaN(nFrames, nPoints, nDims);
    labels = cell(1, nPoints);
    dims   = {'X','Y','Z'};

    for pi = 1:nPoints
        base = pt_bases{pi};
        % Anatomical label if known
        row = strcmp(POINT_ANATOMY(:,1), base);
        if any(row)
            labels{pi} = [base '_' POINT_ANATOMY{row, 2}];
        else
            labels{pi} = base;
        end

        for di = 1:nDims
            candidates = { [base '_' dims{di}], [base '_' lower(dims{di})] };
            for cv = 1:numel(candidates)
                if ismember(candidates{cv}, colNames)
                    pts(:, pi, di) = T.(candidates{cv});
                    break;
                end
            end
        end
    end

    [~, stem] = fileparts(filename);
    fish_points = struct( ...
        'name',        stem, ...
        'frames',      (1:nFrames)', ...
        'point_names', {labels}, ...
        'points',      pts, ...
        'has_z',       has_z, ...
        'format',      'numbered');

    nDimStr = sprintf('%dD', nDims);
    fprintf('Loaded: %s  [%s, %d frames, %d points]  [numbered pt format]\n', ...
            stem, nDimStr, nFrames, nPoints);
    fprintf('  Points: %s\n', strjoin(labels, ' | '));
end
