function fish_points = load_fish_curves(filename)
% LOAD_FISH_CURVES  Load a CURVES-format midline CSV (pre-transformed body coordinates).
%
%   fish_points = load_fish_curves(filename)
%
%   FORMAT E (CURVES):
%     Row 0:  body station positions in mm from snout, every 10 mm (60, 70, ..., 240).
%             Columns come in pairs: [station_mm, NaN, station_mm, NaN, ...]
%     Rows 1+: per-frame midline data.
%             Even columns (0-indexed): X = normalized body-axis position (0=head, ~1=tail, in BL)
%             Odd  columns (0-indexed): Y = lateral displacement in body lengths (BL)
%             19 body stations × 2 columns = 38 columns total.
%
%   The data is ALREADY in body-frame coordinates (equivalent to the output of
%   transform_fish), so no further coordinate transformation is needed.
%   Pass the result directly to compute_kinematics — skip transform_fish.
%
%   Typical usage:
%     fp   = load_fish_curves('..._CURVES_fish1.csv');
%     kine = compute_kinematics(fp, fps, min_freq);
%
%   OUTPUT  fish_points — struct (same schema as transform_fish output):
%     .name        string  — filename stem
%     .frames      [nFrames x 1]
%     .point_names {1 x nPoints}  body station labels, e.g. {'s60mm','s70mm',...}
%     .points      [nFrames x nPoints x 2]  raw (X, Y) columns from file
%     .has_z       false  (2-D only)
%     .format      'curves'
%     .body_pos_mm [1 x nPoints]  body station positions in mm (60–240)
%     .X           [nFrames x nPoints]  normalized body-axis position (pre-computed)
%     .Y           [nFrames x nPoints]  lateral displacement in BL (pre-computed)
%     .pre_transformed  true  — signals that transform_fish must be skipped

    %% ---- Read raw file (no header parsing by readtable — first row is our header) ----
    opts = detectImportOptions(filename);
    opts.VariableNamingRule = 'preserve';
    opts.DataLines          = [1 Inf];   % read everything
    opts.VariableNamesLine  = 0;         % no variable-name line
    opts.EmptyLineRule      = 'read';
    T = readtable(filename, opts);

    raw = table2array(T);   % [nRows x nCols],  row 1 = header, rows 2+ = data

    %% ---- Parse header row ----
    header_row = raw(1, :);
    % Even columns (1-indexed odd): body position in mm; odd columns: NaN
    body_pos_mm = header_row(1:2:end);       % [1 x nStations]
    body_pos_mm = body_pos_mm(~isnan(body_pos_mm));
    nStations   = numel(body_pos_mm);

    %% ---- Extract data rows ----
    data   = raw(2:end, :);       % [nFrames x (2*nStations)]
    nFrames = size(data, 1);

    % Drop trailing all-NaN rows (some files have an extra blank row)
    valid_rows = ~all(isnan(data), 2);
    data       = data(valid_rows, :);
    nFrames    = size(data, 1);

    %% ---- Split into X and Y matrices ----
    % Even 1-indexed columns = X (lateral displacement in BL)
    % Odd  1-indexed columns = Y (dorso-ventral or lateral — see notes)
    %
    % NOTE: In this format X is the normalized position along the body axis
    % (increases 0→1 head to tail), and Y is lateral displacement in BL.
    % Column ordering: col1=X@s1, col2=Y@s1, col3=X@s2, col4=Y@s2, ...
    X_mat = data(:, 1:2:end);   % [nFrames x nStations]
    Y_mat = data(:, 2:2:end);

    % Guard: if column count doesn't match body_pos_mm, trim
    n_avail = min([size(X_mat,2), size(Y_mat,2), nStations]);
    X_mat       = X_mat(:, 1:n_avail);
    Y_mat       = Y_mat(:, 1:n_avail);
    body_pos_mm = body_pos_mm(1:n_avail);
    nStations   = n_avail;

    %% ---- Build point names ----
    labels = arrayfun(@(mm) sprintf('s%dmm', mm), body_pos_mm, 'UniformOutput', false);

    %% ---- Pack points array (for compatibility with pipeline) ----
    pts       = NaN(nFrames, nStations, 2);
    pts(:,:,1) = X_mat;
    pts(:,:,2) = Y_mat;

    %% ---- Output ----
    [~, stem] = fileparts(filename);

    fish_points = struct( ...
        'name',             stem, ...
        'frames',           (1:nFrames)', ...
        'point_names',      {labels}, ...
        'points',           pts, ...
        'has_z',            false, ...
        'format',           'curves', ...
        'body_pos_mm',      body_pos_mm, ...
        'X',                X_mat, ...
        'Y',                Y_mat, ...
        'pre_transformed',  true);

    fprintf('Loaded CURVES: %s  [%d frames, %d body stations, %.0f–%.0f mm]\n', ...
            stem, nFrames, nStations, body_pos_mm(1), body_pos_mm(end));
    fprintf('  X range: [%.3f, %.3f] BL   Y range: [%.3f, %.3f] BL\n', ...
            min(X_mat(:)), max(X_mat(:)), min(Y_mat(:)), max(Y_mat(:)));
end
