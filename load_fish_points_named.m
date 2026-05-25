function fish_points = load_fish_points_named(filename, selected_points, point_order)
% LOAD_FISH_POINTS_NAMED  Load named-column or dual-camera 2D tracking CSV.
%
%   fish_points = load_fish_points_named(filename, selected_points, point_order)
%   fish_points = load_fish_points_named(filename)   % prints available points and returns
%
%   Supported column formats:
%     FORMAT B (named):       eye_X, eye_Y[, eye_Z], snout_X, ...
%     FORMAT D (dual-camera): pt1_cam1_X, pt1_cam1_Y, pt1_cam2_X, pt1_cam2_Y, ...
%
%   FORMAT D anatomy:
%     pt1  = base of pectoral fin (right)     pt7  = tip of pelvic fin
%     pt2  = tip of pectoral fin (right)      pt8  = base of anal fin
%     pt3  = peduncle                          pt9  = tip of anal fin
%     pt4  = tip of caudal fin                pt10 = tip of dorsal fin
%     pt5  = eye (right)                      pt11 = eye (left / #2)
%     pt6  = base of pelvic fin               pt12 = pectoral fin #2 (left base)
%     cam1 = lateral view   cam2 = ventral view
%
%   For FORMAT D, the pectoral fin phase analysis (pt2 vs pt12) is the
%   primary output: see .pect_phase_result field.
%
%   INPUTS (FORMAT B):
%     selected_points  - cell array of point base-names to include
%     point_order      - integer ordering vector (head=1 to tail=end)
%
%   OUTPUT  fish_points — 1-element struct:
%     .name          string  — filename stem
%     .frames        [nFrames x 1]
%     .point_names   {1 x nPoints}  ordered labels
%     .points        [nFrames x nPoints x nDims]   nDims = 2 or 3
%     .has_z         logical
%     .format        string: 'named' or 'dual_camera'
%     .cam_data      struct  (FORMAT D only — per-camera per-point raw data)
%     .pect_phase_result  struct  (FORMAT D only — pectoral fin phase)

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

    %% If no selection given — print available names and return
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
                        pair_data(:, ei, di) = T.(colNames{col});
                    end
                end
            end
            pts(:, pi, :) = mean(pair_data, 2, 'omitnan');
            labels{pi}    = strjoin(entry, '+');
        else
            for di = 1:nDims
                col = find_col(colNames, entry, dims{di});
                if ~isempty(col)
                    pts(:, pi, di) = T.(colNames{col});
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
%  FORMAT D LOADER: pt<N>_cam<M>_X / _Y
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

    % Find all pt numbers and camera numbers
    tok_all  = regexp(colNames, '^(pt\d+)_(cam\d+)_[XxYy]$', 'tokens');
    matched  = tok_all(~cellfun(@isempty, tok_all));
    pt_keys  = cellfun(@(t) t{1}{1}, matched, 'UniformOutput', false);
    cam_keys = cellfun(@(t) t{1}{2}, matched, 'UniformOutput', false);

    pt_bases  = unique(pt_keys,  'stable');
    cam_bases = unique(cam_keys, 'stable');

    % Sort numerically
    pt_nums  = cellfun(@(s) str2double(regexp(s,'\d+','match','once')), pt_bases);
    [~, ord] = sort(pt_nums);
    pt_bases = pt_bases(ord);

    cam_nums  = cellfun(@(s) str2double(regexp(s,'\d+','match','once')), cam_bases);
    [~, ord2] = sort(cam_nums);
    cam_bases = cam_bases(ord2);

    nPts  = numel(pt_bases);
    nCams = numel(cam_bases);

    % Store per-camera per-point data (X, Y) as a struct
    cam_data = struct();
    for ci = 1:nCams
        cam = cam_bases{ci};
        for pi = 1:nPts
            pt = pt_bases{pi};
            col_x = find_col_dual(colNames, pt, cam, 'X');
            col_y = find_col_dual(colNames, pt, cam, 'Y');
            if ~isempty(col_x)
                cam_data.(cam).(pt).X = T.(colNames{col_x});
            else
                cam_data.(cam).(pt).X = NaN(nFrames, 1);
            end
            if ~isempty(col_y)
                cam_data.(cam).(pt).Y = T.(colNames{col_y});
            else
                cam_data.(cam).(pt).Y = NaN(nFrames, 1);
            end
        end
    end

    % Build combined points array for kinematics:
    % Use cam1 (lateral) as primary — take X/Y for each point
    % Points = those that have labels matching anatomy
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
        % Use cam1 (lateral view) as default for kinematics
        if nCams >= 1
            cam = cam_bases{1};
            pts(:, pi, 1) = cam_data.(cam).(pt).X;
            pts(:, pi, 2) = cam_data.(cam).(pt).Y;
        end
    end

    % ---- Pectoral fin phase analysis (pt2 vs pt12) ----
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
%  PECTORAL FIN PHASE ANALYSIS (FORMAT D)
% =========================================================================
function result = compute_pect_phase(cam_data, cam_bases, pt_bases, nFrames)
% Computes the phase relationship between pectoral fin tip right (pt2)
% and pectoral fin base left (pt12) using both camera views.
%
% Uses the cam2 (ventral) view for the most sensitive dorso-ventral signal,
% falls back to cam1 if cam2 is unavailable for a given point.
%
% Classification:
%   In-phase:        |phase_shift| <= 45 deg or >= 315 deg
%   Antiphase:       135 <= |phase_shift| <= 225 deg
%   Intermediate:    all other values

    result = struct();

    % Prefer cam2 (ventral) for pectoral fin phase (better dorsoventral resolution)
    % Fall back to cam1 if cam2 not available
    function [sig, cam_used] = best_signal(pt)
        cam_used = '';
        sig = [];
        % Try cam2 first (ventral)
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

    % Use only frames where both are valid
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

    % FFT cross-correlation to find dominant frequency and phase
    N = n_valid;
    F2  = fft(s2);
    F12 = fft(s12);
    pow2  = abs(F2(1:floor(N/2)+1)).^2;
    freqs = (0:floor(N/2)) / N;   % normalized frequency

    % Find dominant beat frequency (ignore DC, i.e., freq > 0.01)
    valid_f = freqs > 0.01;
    [~, fi] = max(pow2(valid_f));
    dom_idx = find(valid_f, 1) + fi - 1;  % index in full FFT

    % Phase of each signal at the dominant frequency
    phase2  = angle(F2(dom_idx))  * 180 / pi;
    phase12 = angle(F12(dom_idx)) * 180 / pi;

    phase_shift = mod(phase12 - phase2, 360);

    % Classification
    if phase_shift <= 45 || phase_shift >= 315
        classification = 'In-phase';
    elseif phase_shift >= 135 && phase_shift <= 225
        classification = 'Antiphase';
    else
        classification = 'Intermediate';
    end

    % Also compute cross-correlation for the time-domain view
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
% Simple normalized cross-correlation (no toolbox required)
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
