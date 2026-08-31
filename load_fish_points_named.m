function fish_points = load_fish_points_named(filename, selected_points, point_order)
% LOAD_FISH_POINTS_NAMED  Load named-column or dual-camera 2D tracking CSV.
%
%   (See original docstring, unchanged. This version patches the FORMAT B
%   point-extraction loop to handle columns MATLAB imported as text
%   because every value in them was "NA", a fully-occluded landmark.
%   Previously this would crash with "Conversion to double from cell is
%   not possible"; now it becomes a clean all-NaN column with a warning.)

    %% Read file
    opts = detectImportOptions(filename);
    opts.VariableNamingRule = 'preserve';
    T = readtable(filename, opts);
    colNames = T.Properties.VariableNames;
    nFrames  = height(T);
    frames   = (1:nFrames)';

    %% ---- Detect FORMAT D: pt<N>_cam<M>_X ----
    is_dualcam = any(~cellfun(@isempty, regexp(colNames, '^pt\d+_cam\d+_[XxYy]$')));

    if is_dualcam
        fish_points = load_dual_camera(T, colNames, nFrames, frames, filename);
        return;
    end

    %% ---- FORMAT B: named columns ----
    tok_x = regexp(colNames, '^(.+)_[Xx]$', 'tokens');
    tok_z = regexp(colNames, '^(.+)_[Zz]$', 'tokens');
    has_z = any(~cellfun(@isempty, tok_z));

    base_names = cellfun(@(t) t{1}{1}, tok_x(~cellfun(@isempty,tok_x)), 'UniformOutput', false);
    base_names = unique(base_names, 'stable');

    %% If no selection given, print available names and return
    if nargin < 2 || isempty(selected_points)
        fprintf('\nAvailable points in %s:\n', filename);
        for i = 1:numel(base_names)
            fprintf('  %2d.  %s\n', i, base_names{i});
        end
        fprintf('\nHas Z dimension: %s\n\n', mat2str(has_z));
        fish_points = struct('name', filename, 'frames', frames, ...
                             'point_names', {base_names}, 'points', [], 'has_z', has_z, ...
                             'format', 'named');
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
            pair_data = NaN(nFrames, numel(entry), nDims);
            for ei = 1:numel(entry)
                for di = 1:nDims
                    col = find_col(colNames, entry{ei}, dims{di});
                    if ~isempty(col)
                        pair_data(:, ei, di) = to_numeric_col(T.(colNames{col}), entry{ei}, dims{di});
                    end
                end
            end
            pts(:, pi, :) = mean(pair_data, 2, 'omitnan');
            labels{pi}    = strjoin(entry, '+');
        else
            for di = 1:nDims
                col = find_col(colNames, entry, dims{di});
                if ~isempty(col)
                    pts(:, pi, di) = to_numeric_col(T.(colNames{col}), entry, dims{di});
                end
            end
            labels{pi} = entry;
        end
    end

    [~, stem] = fileparts(filename);
    fish_points = struct( ...
        'name',        stem, ...
        'frames',      frames, ...
        'point_names', {labels}, ...
        'points',      pts, ...
        'has_z',       has_z, ...
        'format',      'named');

    nDimStr = sprintf('%dD', nDims);
    fprintf('Loaded: %s  [%s, %d frames, %d points]  [named format]\n', stem, nDimStr, nFrames, nPoints);
    fprintf('  Points: %s\n', strjoin(labels, ' → '));
end


% =========================================================================
%  NEW HELPER: coerce a possibly-text column (from a 100%-NA landmark)
%  to numeric, warning clearly instead of letting a crash happen later.
% =========================================================================
function v = to_numeric_col(col_data, pointName, dimLabel)
    if iscell(col_data) || isstring(col_data)
        col_data = str2double(col_data);
        if all(isnan(col_data))
            warning(['load_fish_points_named: "%s_%s" is 100%% missing in this file ' ...
                     '(fully occluded/untracked landmark), likely the far side of the ' ...
                     'animal from the camera. Consider using the corresponding left/near-' ...
                     'side point instead.'], pointName, dimLabel);
        end
    end
    v = col_data;
end


% =========================================================================
%  FORMAT D LOADER: pt<N>_cam<M>_X / _Y   (unchanged from original except
%  for the same to_numeric_col guard on cam_data assignments)
% =========================================================================
function fish_points = load_dual_camera(T, colNames, nFrames, frames, filename)

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

    tok_all  = regexp(colNames, '^(pt\d+)_(cam\d+)_[XxYy]$', 'tokens');
    matched  = tok_all(~cellfun(@isempty, tok_all));
    pt_keys  = cellfun(@(t) t{1}{1}, matched, 'UniformOutput', false);
    cam_keys = cellfun(@(t) t{1}{2}, matched, 'UniformOutput', false);

    pt_bases  = unique(pt_keys,  'stable');
    cam_bases = unique(cam_keys, 'stable');

    pt_nums  = cellfun(@(s) str2double(regexp(s,'\d+','match','once')), pt_bases);
    [~, ord] = sort(pt_nums);
    pt_bases = pt_bases(ord);

    cam_nums  = cellfun(@(s) str2double(regexp(s,'\d+','match','once')), cam_bases);
    [~, ord2] = sort(cam_nums);
    cam_bases = cam_bases(ord2);

    nPts  = numel(pt_bases);
    nCams = numel(cam_bases);

    cam_data = struct();
    for ci = 1:nCams
        cam = cam_bases{ci};
        for pi = 1:nPts
            pt = pt_bases{pi};
            col_x = find_col_dual(colNames, pt, cam, 'X');
            col_y = find_col_dual(colNames, pt, cam, 'Y');
            if ~isempty(col_x)
                cam_data.(cam).(pt).X = to_numeric_col(T.(colNames{col_x}), [pt '_' cam], 'X');
            else
                cam_data.(cam).(pt).X = NaN(nFrames, 1);
            end
            if ~isempty(col_y)
                cam_data.(cam).(pt).Y = to_numeric_col(T.(colNames{col_y}), [pt '_' cam], 'Y');
            else
                cam_data.(cam).(pt).Y = NaN(nFrames, 1);
            end
        end
    end

    nDims  = 2;
    pts    = NaN(nFrames, nPts, nDims);
    labels = cell(1, nPts);

    for pi = 1:nPts
        pt = pt_bases{pi};
        row = strcmp(POINT_ANATOMY(:,1), pt);
        if any(row)
            labels{pi} = [pt '_' POINT_ANATOMY{row, 2}];
        else
            labels{pi} = pt;
        end
        if nCams >= 1
            cam = cam_bases{1};
            pts(:, pi, 1) = cam_data.(cam).(pt).X;
            pts(:, pi, 2) = cam_data.(cam).(pt).Y;
        end
    end

    pect_phase = compute_pect_phase(cam_data, cam_bases, pt_bases, nFrames);

    [~, stem] = fileparts(filename);
    fish_points = struct( ...
        'name',               stem, ...
        'frames',             frames, ...
        'point_names',        {labels}, ...
        'points',             pts, ...
        'has_z',              false, ...
        'format',             'dual_camera', ...
        'cam_data',           cam_data, ...
        'cam_names',          {cam_bases}, ...
        'pt_bases',           {pt_bases}, ...
        'pect_phase_result',  pect_phase);

    fprintf('Loaded: %s  [2D dual-camera, %d frames, %d points, %d cameras]  [dual-camera format]\n', ...
            stem, nFrames, nPts, nCams);
    fprintf('  Cameras: %s\n', strjoin(cam_bases, ', '));
    fprintf('  Points:  %s\n', strjoin(pt_bases,  ' | '));
    if ~isempty(pect_phase)
        fprintf('  Pectoral phase: %s (phase shift = %.1f deg)\n', ...
                pect_phase.classification, pect_phase.phase_shift_deg);
    end
end


% =========================================================================
%  PECTORAL FIN PHASE ANALYSIS (FORMAT D), unchanged from original
% =========================================================================
function result = compute_pect_phase(cam_data, cam_bases, pt_bases, nFrames)
    result = struct();

    function [sig, cam_used] = best_signal(pt)
        cam_used = '';
        sig = [];
        for ci = numel(cam_bases):-1:1
            cam = cam_bases{ci};
            if isfield(cam_data, cam) && isfield(cam_data.(cam), pt)
                y = cam_data.(cam).(pt).Y;
                if sum(~isnan(y)) > nFrames * 0.3
                    sig = y;
                    cam_used = cam;
                    return;
                end
            end
        end
    end

    has_pt2  = ismember('pt2',  pt_bases);
    has_pt12 = ismember('pt12', pt_bases);

    if ~has_pt2 || ~has_pt12
        result.classification = 'N/A (pt2 or pt12 missing)';
        result.phase_shift_deg = NaN;
        result.n_valid = 0;
        return;
    end

    [sig2,  cam2_used]  = best_signal('pt2');
    [sig12, cam12_used] = best_signal('pt12');

    if isempty(sig2) || isempty(sig12)
        result.classification = 'N/A (insufficient data)';
        result.phase_shift_deg = NaN;
        result.n_valid = 0;
        return;
    end

    valid = ~isnan(sig2) & ~isnan(sig12);
    n_valid = sum(valid);

    if n_valid < 20
        result.classification = 'N/A (too few valid frames)';
        result.phase_shift_deg = NaN;
        result.n_valid = n_valid;
        return;
    end

    s2  = sig2(valid)  - mean(sig2(valid));
    s12 = sig12(valid) - mean(sig12(valid));

    N = n_valid;
    F2  = fft(s2);
    F12 = fft(s12);
    pow2  = abs(F2(1:floor(N/2)+1)).^2;
    freqs = (0:floor(N/2)) / N;

    valid_f = freqs > 0.01;
    [~, fi] = max(pow2(valid_f));
    dom_idx = find(valid_f, 1) + fi - 1;

    phase2  = angle(F2(dom_idx))  * 180 / pi;
    phase12 = angle(F12(dom_idx)) * 180 / pi;

    phase_shift = mod(phase12 - phase2, 360);

    if phase_shift <= 45 || phase_shift >= 315
        classification = 'In-phase';
    elseif phase_shift >= 135 && phase_shift <= 225
        classification = 'Antiphase';
    else
        classification = 'Intermediate';
    end

    [xcorr_vals, lags] = xcorr_simple(s2, s12);
    [~, max_lag_idx]   = max(xcorr_vals);
    peak_lag           = lags(max_lag_idx);

    result.classification  = classification;
    result.phase_shift_deg = phase_shift;
    result.phase_right_deg = phase2;
    result.phase_left_deg  = phase12;
    result.dom_freq_norm   = freqs(dom_idx);
    result.peak_lag_frames = peak_lag;
    result.n_valid         = n_valid;
    result.sig_right       = sig2;
    result.sig_left        = sig12;
    result.valid_mask      = valid;
    result.xcorr_vals      = xcorr_vals;
    result.xcorr_lags      = lags;
    result.cam_right       = cam2_used;
    result.cam_left        = cam12_used;
end


function [c, lags] = xcorr_simple(x, y)
    N    = length(x);
    maxl = min(N-1, round(N/2));
    lags = -maxl:maxl;
    c    = zeros(1, numel(lags));
    sx   = std(x, 'omitnan');
    sy   = std(y, 'omitnan');
    denom = (N * sx * sy);
    if denom == 0, return; end
    for k = 1:numel(lags)
        lag = lags(k);
        if lag >= 0
            xi = x(1:N-lag);
            yi = y(lag+1:N);
        else
            xi = x(-lag+1:N);
            yi = y(1:N+lag);
        end
        c(k) = sum((xi - mean(xi)) .* (yi - mean(yi)), 'omitnan') / denom;
    end
end


% -------------------------------------------------------------------------
function idx = find_col(colNames, base, dim)
    pattern = sprintf('^%s_%s$', regexptranslate('escape', base), dim);
    idx = find(~cellfun(@isempty, regexpi(colNames, pattern)), 1);
end

function idx = find_col_dual(colNames, pt, cam, dim)
    pattern = sprintf('^%s_%s_%s$', regexptranslate('escape', pt), ...
                      regexptranslate('escape', cam), dim);
    idx = find(~cellfun(@isempty, regexpi(colNames, pattern)), 1);
end
