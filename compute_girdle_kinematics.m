function girdle = compute_girdle_kinematics(csvPath, girdlePointName, fish_points, fps, min_freq)
% COMPUTE_GIRDLE_KINEMATICS  Protraction/retraction kinematics of a
%                             pectoral or pelvic girdle base point.
%
%   girdle = compute_girdle_kinematics(csvPath, girdlePointName, fish_points, fps, min_freq)
%
%   Expresses the girdle base point in the SAME body-relative coordinate
%   frame the midline uses (via apply_body_transform + fish_points'
%   .transform_params), then reports how much the girdle itself moves
%   fore-aft (protraction/retraction) and side-to-side over the swim
%   cycle. This is a genuinely separate signal from the fin's own
%   yaw/pitch (compute_fin_kinematics) — it captures girdle ROTATION
%   about the body itself, not fin sweep relative to the girdle.
%
%   INPUTS
%     csvPath          - CSV containing the girdle point's raw X/Y[/Z]
%                         columns (same file used for load_fish_points_named).
%     girdlePointName   - base column name, e.g. 'RPectBase' (expects
%                          columns girdlePointName_X, _Y[, _Z]).
%     fish_points        - the SAME transform_fish() output used for the
%                          midline (must have .transform_params — requires
%                          the updated transform_fish.m).
%     fps, min_freq       - as elsewhere.
%
%   OUTPUT  girdle — struct:
%     .X / .Y / .Z            [nFrames x 1]  girdle position in body-relative
%                               BL-normalized coordinates (X: fore-aft
%                               position along the body axis — protraction
%                               moves this forward/smaller, retraction
%                               moves it aft/larger; Y: lateral splay)
%     .protraction_range_BL    fore-aft (X) excursion range
%     .lateral_range_BL        side-to-side (Y) excursion range
%     .girdle_freq_Hz          dominant oscillation frequency of the
%                               fore-aft (X) signal — the girdle's own
%                               protraction/retraction beat rate
%     .n_valid / .n_frames / .pct_valid
%
%   REQUIRES: apply_body_transform.m on the path.

    if numel(fish_points) ~= 1
        error('compute_girdle_kinematics: pass a single-animal fish_points element, e.g. fish_points(1).');
    end
    if ~isfield(fish_points, 'transform_params')
        error(['compute_girdle_kinematics: fish_points.transform_params not found. ' ...
               'Re-run transform_fish() with the updated version that exports ' ...
               'transform_params before calling this function.']);
    end

    %% ---- Load the raw girdle point from the CSV ----
    opts = detectImportOptions(csvPath);
    opts.VariableNamingRule = 'preserve';
    T = readtable(csvPath, opts);
    cols = T.Properties.VariableNames;
    nFrames = height(T);

    has_z = ismember([girdlePointName '_Z'], cols) || ismember([girdlePointName '_z'], cols);
    raw_xyz = NaN(nFrames, 2 + has_z);

    dim_labels = {'X','Y','Z'};
    for d = 1:(2+has_z)
        cand = {[girdlePointName '_' dim_labels{d}], [girdlePointName '_' lower(dim_labels{d})]};
        found = false;
        for c = 1:numel(cand)
            if ismember(cand{c}, cols)
                col_data = T.(cand{c});

                % CHANGE NOTE (bug fix): if a landmark is 100% missing in
                % this file (e.g. the far-side fin, occluded from camera —
                % same pattern as BP_16_RPectoral1 etc.), MATLAB's readtable
                % can't infer the column is numeric and imports it as text
                % (cell/char) instead of double. Convert explicitly —
                % str2double() turns unparseable text (including "NA")
                % into NaN, which is what we want: a fully-occluded point
                % becomes a clean all-NaN column, not a crash.
                if iscell(col_data) || isstring(col_data)
                    col_data = str2double(col_data);
                    if all(isnan(col_data))
                        warning(['compute_girdle_kinematics: "%s" is 100%% missing in this ' ...
                                 'file (fully occluded/untracked landmark) — likely the far ' ...
                                 'side of the animal from the camera. Try the corresponding ' ...
                                 'left/near-side point instead.'], cand{c});
                    end
                end

                % Defensive check: a valid MATLAB table column should
                % always have exactly nFrames rows, but if this ever
                % trips, fail with a message that actually says why,
                % instead of MATLAB's generic "different number of
                % elements" assignment error.
                if numel(col_data) ~= nFrames
                    error(['compute_girdle_kinematics: column "%s" has %d values but the ' ...
                           'file has %d rows (class: %s). This points to a malformed/ragged ' ...
                           'CSV row somewhere in the source file — check it opens cleanly ' ...
                           'in Excel/a text editor with a consistent column count on every row.'], ...
                           cand{c}, numel(col_data), nFrames, class(T.(cand{c})));
                end

                raw_xyz(:,d) = col_data(:);
                found = true; break;
            end
        end
        if ~found
            error('compute_girdle_kinematics: column "%s_%s" not found in %s.', ...
                   girdlePointName, dim_labels{d}, csvPath);
        end
    end

    if nFrames ~= numel(fish_points.transform_params)
        error(['compute_girdle_kinematics: frame count mismatch — %s has %d rows but ' ...
               'fish_points.transform_params has %d frames. Make sure this is the same ' ...
               'source file/trial used for the midline.'], csvPath, nFrames, numel(fish_points.transform_params));
    end

    %% ---- Project into body-relative frame ----
    [X, Y, Z] = apply_body_transform(raw_xyz, fish_points.transform_params);

    n_valid  = sum(~isnan(X) & ~isnan(Y));
    pct_valid = 100 * n_valid / nFrames;

    if n_valid == 0
        warning(['compute_girdle_kinematics: %s has ZERO valid frames after projection — ' ...
                 'girdle metrics will be NaN. Check that %s is tracked in at least some ' ...
                 'frames where the midline transform also succeeded.'], girdlePointName, girdlePointName);
        protraction_range = NaN; lateral_range = NaN; girdle_freq = NaN;
    else
        protraction_range = max(X(~isnan(X))) - min(X(~isnan(X)));
        lateral_range      = max(Y(~isnan(Y))) - min(Y(~isnan(Y)));

        % Fore-aft oscillation frequency — same NaN-safe dominant_freq logic
        % as compute_kinematics.m (all-NaN / zero-variance -> NaN, not a
        % fabricated bin).
        X_filled = fill_nan_local(X);
        girdle_freq = dominant_freq_local(X_filled, fps, min_freq);
    end

    girdle.name                = girdlePointName;
    girdle.X = X; girdle.Y = Y; girdle.Z = Z;
    girdle.protraction_range_BL = protraction_range;
    girdle.lateral_range_BL      = lateral_range;
    girdle.girdle_freq_Hz        = girdle_freq;
    girdle.n_valid  = n_valid;
    girdle.n_frames = nFrames;
    girdle.pct_valid = pct_valid;

    fprintf(['Girdle (%s): %d/%d frames (%.1f%%) | protraction range=%s BL  ' ...
             'lateral range=%s BL  freq=%s Hz\n'], girdlePointName, n_valid, nFrames, pct_valid, ...
            fmt_val(protraction_range), fmt_val(lateral_range), fmt_val(girdle_freq));
end


function s = fmt_val(v)
    if isnan(v), s = 'NaN'; else, s = sprintf('%.4f', v); end
end

function y = fill_nan_local(y)
    t = (1:length(y))';
    good = ~isnan(y);
    if sum(good) < 2, y(:) = NaN; return; end
    y(~good) = interp1(t(good), y(good), t(~good), 'linear', 'extrap');
end

function f_dom = dominant_freq_local(y, fs, min_freq)
    if all(isnan(y)) || (max(y) - min(y)) < eps
        f_dom = NaN; return;
    end
    N = length(y);
    Y = fft(y - mean(y));
    power = (2/N) * abs(Y(1:floor(N/2)+1)).^2;
    freqs = fs * (0:floor(N/2)) / N;
    valid = freqs >= min_freq;
    if ~any(valid), f_dom = NaN; return; end
    power_valid = power(valid);
    freqs_valid = freqs(valid);
    [~, idx] = max(power_valid);
    f_dom = freqs_valid(idx);
end