function fish_points = load_fish_points(filename)
% LOAD_FISH_POINTS  Load a DLC-style tracking CSV with columns Fish1_P1_x/y/z.
%
%   fish_points = load_fish_points(filename)
%
%   INPUT
%     filename    - path to CSV with columns: frame, FishN_Pk_x, FishN_Pk_y[, FishN_Pk_z]
%
%   OUTPUT  fish_points — struct array, one element per fish:
%     .name        string,  e.g. 'Fish1'
%     .frames      [nFrames x 1]
%     .point_names {1 x nPoints}  cell of point labels, e.g. {'P1','P2',...}
%     .points      [nFrames x nPoints x nDims]   nDims = 2 (xy) or 3 (xyz)
%     .has_z       logical scalar
%
%   All downstream functions (transform_fish, compute_kinematics) accept any
%   number of points and work in 2-D or 3-D automatically.

    opts = detectImportOptions(filename);
    opts.VariableNamingRule = 'preserve';
    T = readtable(filename, opts);

    colNames = T.Properties.VariableNames;
    frames   = T.frame;
    nFrames  = height(T);

    % Match columns: FishN_Pk_x / _y / _z
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

    fish_points(nFish) = struct('name','','frames',[],'point_names',{{}},'points',[],'has_z',false);

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
    end

    nDimStr = sprintf('%dD', nDims);
    fprintf('Loaded: %s  [%s, %d frames, %d fish, %d points]\n', ...
            filename, nDimStr, nFrames, nFish, nPoints);
end
