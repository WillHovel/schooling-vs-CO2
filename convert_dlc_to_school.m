function convert_dlc_to_school(input)
%CONVERT_DLC_TO_SCHOOL  Convert wide DeepLabCut exports to Kinemetrix School format.
%
%   Converts a DeepLabCut "wide" export whose columns are named
%       index, ind<N>_<bodypart>_<coord>, ...
%   into the Kinemetrix School Metrics convention
%       frame, Fish<N>_P<k>_<coord>, ...
%   which load_fish_points() recognises. This is the MATLAB, multi-file
%   replacement for the single-file convert_f0_to_school.py script.
%
%   Everything is discovered from each file's header, so the converter works
%   on ANY file of this format rather than one hard-coded path:
%     * individuals     - the set of N in ind<N>, sorted numerically;
%     * bodyparts       - discovered in order of first appearance, assigned
%                         point numbers P1..Pk in that order;
%     * coordinates     - discovered in order of first appearance (x, y, z).
%
%   'NA' / 'N/A' / empty cells are written as 'NaN' so readtable keeps every
%   coordinate column numeric (the same convention as the Python script).
%
%   USAGE
%     convert_dlc_to_school                  prompts for a folder, converts all *.csv
%     convert_dlc_to_school(DIR)             converts every *.csv in DIR
%     convert_dlc_to_school(FILE)            converts a single file
%     convert_dlc_to_school({F1,F2,...})     converts the listed files
%
%   Each input produces <stem>_school.csv next to the source file.

    if nargin < 1 || isempty(input)
        input = uigetdir(pwd, 'Select the folder of DLC CSV files');
        if isequal(input, 0)
            return;
        end
    end

    files = resolve_files(input);
    if isempty(files)
        warning('convert_dlc_to_school:noFiles', 'No CSV files matched.');
        return;
    end

    for i = 1:numel(files)
        try
            convert_one(files{i});
        catch ME
            warning('convert_dlc_to_school:convertFailed', ...
                'Skipping %s: %s', files{i}, ME.message);
        end
    end
end

% -------------------------------------------------------------------------

function files = resolve_files(input)
    if ischar(input) || isstring(input)
        input = char(input);
        if isfolder(input)
            d = dir(fullfile(input, '*.csv'));
            names = {d.name};
            names = names(~endsWith(names, '_school.csv'));   % skip prior outputs
            files = fullfile(input, names);
        elseif isfile(input)
            files = {input};
        else
            files = {};
        end
    elseif iscell(input)
        files = input;
    else
        error('convert_dlc_to_school:badInput', ...
            'Input must be a folder, a file path, or a cell array of file paths.');
    end
    files = files(:)';
end

% -------------------------------------------------------------------------

function convert_one(file)
    [p, n, ~] = fileparts(file);
    outfile = fullfile(p, [n '_school.csv']);

    fid = fopen(file, 'r');
    if fid < 0
        error('convert_dlc_to_school:open', 'Could not open %s', file);
    end
    c = onCleanup(@() fclose(fid));

    % ---- Header ----
    headerLine = fgetl(fid);
    if ~ischar(headerLine) || isempty(strtrim(headerLine))
        error('convert_dlc_to_school:empty', 'Empty or missing header in %s', file);
    end
    headerLine = regexprep(headerLine, ['^\xEF\xBB\xBF'], '');   % strip BOM if present
    header = strsplit(headerLine, ',');
    header = strtrim(header);
    header = regexprep(header, '^"|"$', '');                     % strip surrounding quotes
    ncols = numel(header);

    srcIdx = containers.Map(header, num2cell(1:ncols));

    % ---- Discover individuals / bodyparts / coordinates from the header ----
    tok = regexp(header, '^ind(\d+)_(.+)_([xyz])$', 'tokens', 'once', 'ignorecase');
    isCoord = ~cellfun(@isempty, tok);
    if ~any(isCoord)
        error('convert_dlc_to_school:format', ...
            'No ind<N>_<bodypart>_<coord> columns in %s; not a wide DLC export.', file);
    end

    % frame column: prefer a column literally named "index" or "frame",
    % otherwise fall back to the first non-coordinate column.
    frameCol = [];
    for ci = 1:ncols
        if ~isCoord(ci)
            nm = lower(header{ci});
            if strcmp(nm, 'index') || strcmp(nm, 'frame')
                frameCol = ci;
                break;
            end
        end
    end
    if isempty(frameCol)
        frameCol = find(~isCoord, 1, 'first');
    end
    if isempty(frameCol)
        error('convert_dlc_to_school:format', 'No frame (index) column found in %s.', file);
    end

    indNum = zeros(1, ncols);
    bodyparts = {};
    coords = {};
    for ci = 1:ncols
        if ~isCoord(ci), continue; end
        t = tok{ci};
        indNum(ci) = str2double(t{1});
        bp = t{2};
        cd = lower(t{3});
        if ~any(strcmp(bodyparts, bp))
            bodyparts{end+1} = bp; %#ok<AGROW>
        end
        if ~any(strcmp(coords, cd))
            coords{end+1} = cd; %#ok<AGROW>
        end
    end
    individuals = unique(indNum(indNum ~= 0));
    nInd = numel(individuals);
    nBp  = numel(bodyparts);
    nCrd = numel(coords);

    % ---- Verify every source column exists before doing any work ----
    missing = {};
    for ii = 1:nInd
        for k = 1:nBp
            for cc = 1:nCrd
                srcName = sprintf('ind%d_%s_%s', individuals(ii), bodyparts{k}, coords{cc});
                if ~isKey(srcIdx, srcName)
                    missing{end+1} = srcName; %#ok<AGROW>
                end
            end
        end
    end
    if ~isempty(missing)
        shown = strjoin(missing(1:min(5, numel(missing))), ', ');
        error('convert_dlc_to_school:missing', ...
            'Missing source column(s) in %s: %s ...', file, shown);
    end

    % ---- Data ----
    fmt = repmat('%q', 1, ncols);
    C = textscan(fid, fmt, 'Delimiter', ',', 'CollectOutput', true);
    data = C{1};                                     % nrows x ncols cell of char
    data = cellfun(@normalize_missing, data, 'UniformOutput', false);
    nrows = size(data, 1);

    % ---- Output header ----
    outHeader = cell(1, 1 + nInd*nBp*nCrd);
    outHeader{1} = 'frame';
    oc = 1;
    for ii = 1:nInd
        for k = 1:nBp
            for cc = 1:nCrd
                oc = oc + 1;
                outHeader{oc} = sprintf('Fish%d_P%d_%s', individuals(ii), k, coords{cc});
            end
        end
    end

    % ---- Reorder data into output column order ----
    outData = cell(nrows, numel(outHeader));
    outData(:, 1) = data(:, frameCol);
    oc = 1;
    for ii = 1:nInd
        for k = 1:nBp
            for cc = 1:nCrd
                srcName = sprintf('ind%d_%s_%s', individuals(ii), bodyparts{k}, coords{cc});
                oc = oc + 1;
                outData(:, oc) = data(:, srcIdx(srcName));
            end
        end
    end

    % ---- Write ----
    fout = fopen(outfile, 'w');
    if fout < 0
        error('convert_dlc_to_school:write', 'Could not open %s for writing', outfile);
    end
    c2 = onCleanup(@() fclose(fout));
    fprintf(fout, '%s\n', strjoin(outHeader, ','));
    for r = 1:nrows
        fprintf(fout, '%s\n', strjoin(outData(r, :), ','));
    end

    fprintf('Converted %s\n', file);
    fprintf('  %d rows -> %d columns (%d individuals x %d bodyparts x %d coords)\n', ...
        nrows, numel(outHeader), nInd, nBp, nCrd);
    fprintf('  output: %s\n', outfile);
end

% -------------------------------------------------------------------------

function v = normalize_missing(s)
    t = strtrim(s);
    if isempty(t) || strcmpi(t, 'NA') || strcmpi(t, 'N/A')
        v = 'NaN';
    else
        v = t;
    end
end
