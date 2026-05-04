function FishKinematicsApp()
% FISHKINEMATICSAPP  Interactive UI for fish midline kinematics + fin analysis.
%
%   Supports two CSV formats:
%     FORMAT A (DLC-style): columns Fish1_P1_x, Fish1_P1_y[, Fish1_P1_z] ...
%     FORMAT B (named):     columns eye_X, eye_Y[, eye_Z], snout_X ...
%
%   REQUIRED FILES (must all be on the MATLAB path):
%     load_fish_points.m          — loads DLC-format CSVs
%     load_fish_points_named.m    — loads named-column CSVs
%     transform_fish.m            — rotates midline to body axis
%     compute_kinematics.m        — FFT interpolation + kinematic metrics
%     compute_fin_kinematics.m    — fin vector angles, speed, trajectory
%
%   REQUIRED MATLAB TOOLBOXES:
%     BASE MATLAB ONLY — no additional toolboxes required (R2019b or later).
%
%     Functions used and their source:
%       fft / ifft / real         — base MATLAB (no Signal Processing Toolbox needed)
%       mean(...,'omitnan')       — base MATLAB R2015a+  (replaces nanmean)
%       std(...,'omitnan')        — base MATLAB R2015a+  (replaces nanstd)
%       interp1 / polyfit         — base MATLAB
%       readtable / detectImportOptions — base MATLAB
%       uifigure / uiaxes / uitabgroup / uipanel / uilistbox /
%         uibutton / uilabel / uitextarea / uieditfield /
%         uidropdown / uicheckbox — base MATLAB (App Building, R2016a+)
%       atan2d / range / diff     — base MATLAB
%
%     NOT required:
%       Signal Processing Toolbox  (fft is base MATLAB)
%       Statistics and Machine Learning Toolbox  (nanmean/nanstd replaced)

    %% ---- Ensure all helper .m files are on the path ----
    % Adds the folder containing FishKinematicsApp.m itself, so that
    % compute_fin_kinematics.m, compute_kinematics.m, transform_fish.m, etc.
    % are always found regardless of MATLAB's current working directory.
    appDir = fileparts(mfilename('fullpath'));
    if ~isempty(appDir)
        addpath(appDir);
    end

    %% ---- Shared state ----
    app.fp        = [];
    app.kine      = [];
    app.fmt       = '';
    app.avail_pts = {};
    app.sel_order = {};
    app.fin_data  = [];   % computed fin analysis results

    %% ---- Figure ----
    fig = uifigure('Name',     'Fish Kinematics Analyser', ...
                   'Position', [60 40 1260 820], ...
                   'Color',    [0.95 0.95 0.95]);

    %% ================================================================
    %  LEFT PANEL  (inputs, point selector, results)
    %% ================================================================
    LP = uipanel(fig, 'Position', [10 10 340 800], ...
                 'BackgroundColor', [1 1 1], 'BorderType', 'line', 'Title', '');

    y = 768;   % running y cursor (top -> bottom)

    % ---- File ----
    y = section_label(LP, 'INPUT FILE', y);
    app.fileField = uieditfield(LP, 'text', ...
        'Position', [12 y-26 240 26], 'Placeholder', 'No file selected', 'Editable', 'off');
    uibutton(LP, 'Text', 'Browse...', 'Position', [258 y-26 70 26], ...
        'ButtonPushedFcn', @(~,~) onBrowse());
    y = y - 40;

    % ---- Format indicator ----
    app.fmtLabel = uilabel(LP, 'Text', '', ...
        'Position', [12 y-20 316 20], 'FontSize', 11, 'FontColor', [0.3 0.5 0.8]);
    y = y - 32;

    % ---- Parameters ----
    y = section_label(LP, 'PARAMETERS', y);
    uilabel(LP, 'Text', 'Frames per second', 'Position', [12 y-24 190 20], 'FontSize', 12);
    app.fpsField = uieditfield(LP, 'numeric', 'Position', [250 y-24 78 26], ...
        'Value', 100, 'Limits', [1 1e5]);
    y = y - 30;
    uilabel(LP, 'Text', 'Min frequency (Hz)', 'Position', [12 y-24 190 20], 'FontSize', 12);
    app.minFreqField = uieditfield(LP, 'numeric', 'Position', [250 y-24 78 26], ...
        'Value', 0.5, 'Limits', [0 1000]);
    y = y - 42;

    % ---- Point selector (shown only for Format B) ----
    app.ptPanel = uipanel(LP, 'Position', [8 y-210 324 210], ...
        'BackgroundColor', [0.97 0.97 1], 'BorderType', 'line', 'Title', 'Point Selection');
    app.ptPanel.Visible = 'off';

    uilabel(app.ptPanel, 'Text', 'Available (check to select)', ...
        'Position', [4 162 160 18], 'FontSize', 10, 'FontWeight', 'bold', 'FontColor', [0.4 0.4 0.4]);

    % Scrollable panel to hold checkboxes — populated dynamically in detectFormat
    app.cbScrollPanel = uipanel(app.ptPanel, 'Position', [4 4 130 158], ...
        'BackgroundColor', [0.97 0.97 1], 'BorderType', 'line', 'Scrollable', 'on');
    app.checkboxes = {};   % cell array of uicheckbox handles, filled in detectFormat

    uibutton(app.ptPanel, 'Text', '-> Add',  'Position', [140 140 70 24], 'ButtonPushedFcn', @(~,~) onAddPoint());
    uibutton(app.ptPanel, 'Text', '<-> Avg', 'Position', [140 112 70 24], 'ButtonPushedFcn', @(~,~) onAddAvg());
    uibutton(app.ptPanel, 'Text', 'Up',      'Position', [140 80  70 24], 'ButtonPushedFcn', @(~,~) onMoveUp());
    uibutton(app.ptPanel, 'Text', 'Down',    'Position', [140 56  70 24], 'ButtonPushedFcn', @(~,~) onMoveDown());
    uibutton(app.ptPanel, 'Text', 'Remove',  'Position', [140 28  70 24], 'ButtonPushedFcn', @(~,~) onRemovePoint());
    uibutton(app.ptPanel, 'Text', 'Clear',   'Position', [140 4   70 24], 'ButtonPushedFcn', @(~,~) onClearPoints());

    uilabel(app.ptPanel, 'Text', 'Selected (head -> tail)', ...
        'Position', [216 162 110 18], 'FontSize', 10, 'FontWeight', 'bold', 'FontColor', [0.4 0.4 0.4]);
    app.selList = uilistbox(app.ptPanel, 'Position', [216 4 100 158], ...
        'Multiselect', 'on', 'Items', {});

    y = y - 218;

    % ---- Run button ----
    y = y - 4;
    app.runBtn = uibutton(LP, 'Text', 'Load & Analyse', ...
        'Position', [12 y-36 316 36], 'FontSize', 13, 'FontWeight', 'bold', ...
        'BackgroundColor', [0.18 0.42 0.75], 'FontColor', [1 1 1], ...
        'ButtonPushedFcn', @(~,~) onRun());
    y = y - 44;

    app.statusLabel = uilabel(LP, 'Text', 'Load a CSV file to begin.', ...
        'Position', [12 y-20 316 20], 'FontSize', 11, 'FontColor', [0.4 0.4 0.4], 'WordWrap', 'on');
    y = y - 30;

    % ---- Animal selector ----
    y = section_label(LP, 'SELECT ANIMAL', y);
    app.fishList = uilistbox(LP, 'Position', [12 y-70 316 68], ...
        'Items', {'(run analysis first)'}, 'ValueChangedFcn', @(~,~) onFishSelected());
    y = y - 78;

    % ---- Results ----
    y = section_label(LP, 'KINEMATIC VALUES', y);
    app.resultsArea = uitextarea(LP, 'Position', [12 10 316 y-14], ...
        'Editable', 'off', 'FontSize', 11, 'FontName', 'Courier New', ...
        'BackgroundColor', [0.97 0.97 0.97], 'Value', {''});

    %% ================================================================
    %  RIGHT PANEL  — Tabbed: Kinematics | Fin Analysis
    %% ================================================================
    tg = uitabgroup(fig, 'Position', [360 10 880 800]);

    % ---- Tab 1: Kinematics ----
    tabKine = uitab(tg, 'Title', 'Kinematics');

    % Midline controls bar
    midCtrlPanel = uipanel(tabKine, 'Position', [4 748 866 44], ...
        'BackgroundColor', [0.93 0.95 1], 'BorderType', 'line', 'Title', '');
    uilabel(midCtrlPanel, 'Text', 'Midlines shown:', ...
        'Position', [8 10 110 22], 'FontSize', 11);
    app.midlineCount = uieditfield(midCtrlPanel, 'numeric', ...
        'Position', [122 10 55 24], 'Value', 20, 'Limits', [1 2000], ...
        'ValueChangedFcn', @(~,~) refreshMidlines());
    uilabel(midCtrlPanel, 'Text', 'Colormap:', ...
        'Position', [192 10 72 22], 'FontSize', 11);
    app.cmapDrop = uidropdown(midCtrlPanel, ...
        'Position', [267 10 100 24], ...
        'Items', {'parula','jet','cool','hot','hsv','winter','spring'}, ...
        'Value', 'parula', ...
        'ValueChangedFcn', @(~,~) refreshMidlines());

    % 2x2 axes on kinematics tab
    titles_k = {'Midlines (FFT interpolated)', 'Amplitude envelope', ...
                 'Curvature profile',           'Beat frequency spectra'};
    pos_k    = {[10 390 420 350]; [440 390 420 350]; [10 10 420 350]; [440 10 420 350]};

    app.ax = gobjects(4,1);
    for i = 1:4
        app.ax(i) = uiaxes(tabKine, 'Position', pos_k{i}, ...
            'BackgroundColor', [1 1 1], 'Box', 'on');
        title(app.ax(i), titles_k{i}, 'FontSize', 11);
        app.ax(i).Toolbar.Visible = 'off';
    end
    xlabel(app.ax(1),'X (BL)');                    ylabel(app.ax(1),'Y (BL)');
    xlabel(app.ax(2),'Body pos (0=head,1=tail)');   ylabel(app.ax(2),'Half-amp (BL)');
    xlabel(app.ax(3),'Body pos (0=head,1=tail)');   ylabel(app.ax(3),'Curvature (1/BL)');
    xlabel(app.ax(4),'Frequency (Hz)');             ylabel(app.ax(4),'Power');
    for i = 1:4
        text(app.ax(i), 0.5, 0.5, 'No data', 'Units','normalized', ...
             'HorizontalAlignment','center','FontSize',12,'Color',[0.75 0.75 0.75]);
    end

    % ---- Tab 2: Fin Analysis ----
    tabFin = uitab(tg, 'Title', 'Fin Analysis (3D)');

    % Fin point selector panel
    finSelPanel = uipanel(tabFin, 'Position', [4 720 866 72], ...
        'BackgroundColor', [0.97 1 0.97], 'BorderType', 'line', 'Title', 'Fin Point Selection');
    uilabel(finSelPanel, 'Text', 'Root point:', ...
        'Position', [8 22 80 20], 'FontSize', 11);
    app.finRootDrop = uidropdown(finSelPanel, 'Position', [92 20 160 24], ...
        'Items', {'(load file first)'}, 'Value', '(load file first)');
    uilabel(finSelPanel, 'Text', 'Tip point:', ...
        'Position', [270 22 75 20], 'FontSize', 11);
    app.finTipDrop = uidropdown(finSelPanel, 'Position', [350 20 160 24], ...
        'Items', {'(load file first)'}, 'Value', '(load file first)');
    uibutton(finSelPanel, 'Text', 'Compute Fin', ...
        'Position', [530 16 130 30], 'FontSize', 11, 'FontWeight', 'bold', ...
        'BackgroundColor', [0.1 0.55 0.3], 'FontColor', [1 1 1], ...
        'ButtonPushedFcn', @(~,~) onComputeFin());
    uibutton(finSelPanel, 'Text', 'Animate', ...
        'Position', [674 16 110 30], 'FontSize', 11, ...
        'BackgroundColor', [0.55 0.3 0.1], 'FontColor', [1 1 1], ...
        'ButtonPushedFcn', @(~,~) onAnimateFin());

    % Fin results text
    app.finResultsArea = uitextarea(tabFin, 'Position', [4 580 866 132], ...
        'Editable', 'off', 'FontSize', 11, 'FontName', 'Courier New', ...
        'BackgroundColor', [0.97 0.97 0.97]);

    % Fin axes: 3 plots
    app.finAx = gobjects(3,1);
    fin_titles = {'Pitch / Roll / Yaw over time', ...
                  'Fin tip distance traveled over time', ...
                  'Fin vector trajectory (tip, XY plane)'};
    fin_pos = {[4 300 580 272]; [596 300 274 272]; [596 4 274 272]};
    fin_xl = {'Frame','Frame','X (mm)'};
    fin_yl = {'Angle (deg)','Cumulative distance','Y (mm)'};
    for i = 1:3
        app.finAx(i) = uiaxes(tabFin, 'Position', fin_pos{i}, ...
            'BackgroundColor', [1 1 1], 'Box', 'on');
        title(app.finAx(i), fin_titles{i}, 'FontSize', 10);
        xlabel(app.finAx(i), fin_xl{i});
        ylabel(app.finAx(i), fin_yl{i});
        app.finAx(i).Toolbar.Visible = 'off';
        text(app.finAx(i), 0.5, 0.5, 'No data', 'Units','normalized', ...
             'HorizontalAlignment','center','FontSize',12,'Color',[0.75 0.75 0.75]);
    end

    % Fin 3D vector axis (bottom-left large area)
    app.finAx3D = uiaxes(tabFin, 'Position', [4 4 580 272], ...
        'BackgroundColor', [0.04 0.04 0.08], 'Box', 'on');
    title(app.finAx3D, 'Fin vector trajectory (3D, colored by frame)', 'FontSize', 10, 'Color', [1 1 1]);
    xlabel(app.finAx3D,'X'); ylabel(app.finAx3D,'Y'); zlabel(app.finAx3D,'Z');
    app.finAx3D.Toolbar.Visible = 'off';
    app.finAx3D.XColor = [0.7 0.7 0.7];
    app.finAx3D.YColor = [0.7 0.7 0.7];
    app.finAx3D.ZColor = [0.7 0.7 0.7];
    text(app.finAx3D, 0.5, 0.5, 'No data', 'Units','normalized', ...
         'HorizontalAlignment','center','FontSize',12,'Color',[0.6 0.6 0.6]);

    %% ================================================================
    %  CALLBACKS
    %% ================================================================

    function onBrowse()
        [fn, fp_] = uigetfile('*.csv', 'Select tracking CSV');
        if isequal(fn, 0), return; end
        fullpath = fullfile(fp_, fn);
        app.fileField.Value = fullpath;
        detectFormat(fullpath);
    end

    function detectFormat(path)
        opts = detectImportOptions(path);
        opts.VariableNamingRule = 'preserve';
        T = readtable(path, opts);
        cols = T.Properties.VariableNames;

        dlc_match = ~cellfun(@isempty, regexp(cols, '^Fish\d+_P\d+_[xyXY]'));
        if any(dlc_match)
            app.fmt = 'DLC';
            app.fmtLabel.Text = 'Format: DLC  (Fish1_P1_x columns)';
            app.ptPanel.Visible = 'off';
            % Fin dropdowns not useful for DLC format
            app.finRootDrop.Items = {'N/A — use named format'};
            app.finTipDrop.Items  = {'N/A — use named format'};
        else
            app.fmt = 'named';
            app.fmtLabel.Text = 'Format: Named  (eye_X, snout_Y ... columns)';
            app.ptPanel.Visible = 'on';
            tok = regexp(cols, '^(.+)_[Xx]$', 'tokens');
            bases = cellfun(@(t) t{1}{1}, tok(~cellfun(@isempty,tok)), 'UniformOutput', false);
            app.avail_pts = bases;
            app.sel_order = {};
            app.selList.Items = {};
            buildCheckboxes(bases);
            % Populate fin dropdowns
            app.finRootDrop.Items = bases;
            app.finTipDrop.Items  = bases;
            if numel(bases) >= 2
                app.finRootDrop.Value = bases{1};
                app.finTipDrop.Value  = bases{2};
            end
        end
    end

    % ---- Build checkbox list in the scroll panel ----
    function buildCheckboxes(bases)
        % Delete any existing checkboxes
        delete(app.cbScrollPanel.Children);
        app.checkboxes = {};
        nPts    = numel(bases);
        rowH    = 22;
        padY    = 4;
        totalH  = nPts * rowH + padY;
        % Make inner panel tall enough to scroll
        app.cbScrollPanel.Position(4) = max(158, totalH);   % outer panel clips to 158
        for i = 1:nPts
            y_cb = totalH - i*rowH;   % top-to-bottom order
            cb = uicheckbox(app.cbScrollPanel, ...
                'Text',     bases{i}, ...
                'Value',    0, ...
                'Position', [4 y_cb 118 rowH-2], ...
                'FontSize', 10);
            app.checkboxes{i} = cb;
        end
        % Scroll panel needs the inner content to be taller than the panel frame;
        % MATLAB auto-enables scrolling when children exceed bounds.
        % Reset the scroll panel height to the clipped size after populating.
        app.cbScrollPanel.Position(4) = 158;
    end

    % ---- Return names of currently checked points ----
    function sel = getCheckedPoints()
        sel = {};
        for i = 1:numel(app.checkboxes)
            if app.checkboxes{i}.Value
                sel{end+1} = app.checkboxes{i}.Text; %#ok<AGROW>
            end
        end
    end

    % ---- Uncheck all checkboxes ----
    function uncheckAll()
        for i = 1:numel(app.checkboxes)
            app.checkboxes{i}.Value = 0;
        end
    end

    % ---- Point list helpers ----
    function onAddPoint()
        sel = getCheckedPoints();
        if isempty(sel)
            uialert(fig, 'Check at least one point to add.', 'No selection'); return
        end
        for k = 1:numel(sel)
            if ~any(strcmp(entry_label(sel{k}), app.selList.Items))
                app.sel_order{end+1} = sel{k};
            end
        end
        uncheckAll();
        refreshSelList();
    end

    function onAddAvg()
        sel = getCheckedPoints();
        if numel(sel) < 2
            uialert(fig, 'Check 2 or more points to average.', 'Average'); return
        end
        lbl = strjoin(sel, '+');
        if ~any(strcmp(lbl, app.selList.Items))
            app.sel_order{end+1} = sel;   % store as cell array of names
        end
        uncheckAll();
        refreshSelList();
    end

    function onMoveUp()
        idx = getSelListIdx();
        if isempty(idx) || idx == 1, return; end
        app.sel_order([idx-1 idx]) = app.sel_order([idx idx-1]);
        refreshSelList();
        app.selList.Value = app.selList.Items{idx-1};
    end

    function onMoveDown()
        idx = getSelListIdx();
        if isempty(idx) || idx == numel(app.sel_order), return; end
        app.sel_order([idx idx+1]) = app.sel_order([idx+1 idx]);
        refreshSelList();
        app.selList.Value = app.selList.Items{idx+1};
    end

    function onRemovePoint()
        idx = getSelListIdx();
        if isempty(idx), return; end
        app.sel_order(idx) = [];
        refreshSelList();
    end

    function onClearPoints()
        app.sel_order = {};
        refreshSelList();
    end

    function refreshSelList()
        lbls = cellfun(@entry_label, app.sel_order, 'UniformOutput', false);
        app.selList.Items = lbls;
    end

    function lbl = entry_label(e)
        if iscell(e), lbl = strjoin(e, '+'); else, lbl = e; end
    end

    function idx = getSelListIdx()
        val = app.selList.Value;
        idx = find(strcmp(val, app.selList.Items), 1);
    end

    % ---- Run ----
    function onRun()
        csvPath = strtrim(app.fileField.Value);
        if isempty(csvPath) || ~isfile(csvPath)
            uialert(fig, 'Select a valid CSV file first.', 'No file'); return
        end
        fps     = app.fpsField.Value;
        minFreq = app.minFreqField.Value;

        setStatus('Loading...');
        try
            if strcmp(app.fmt, 'DLC')
                app.fp = load_fish_points(csvPath);
            else
                if numel(app.sel_order) < 3
                    uialert(fig, 'Select at least 3 points (head, middle(s), tail).', 'Points'); return
                end
                app.fp = load_fish_points_named(csvPath, app.sel_order, []);
                app.fp = app.fp(1);
            end
        catch ME; uialert(fig, ME.message, 'Load error'); setStatus('Load failed.'); return; end

        setStatus('Transforming...');
        try
            app.fp = transform_fish(app.fp);
        catch ME; uialert(fig, ME.message, 'Transform error'); setStatus('Failed.'); return; end

        setStatus('Computing kinematics...');
        try
            app.kine = compute_kinematics(app.fp, fps, minFreq);
        catch ME; uialert(fig, ME.message, 'Kinematics error'); setStatus('Failed.'); return; end

        names = {app.kine.name};
        app.fishList.Items = names;
        app.fishList.Value = names{1};
        setStatus(sprintf('Done. %d animal(s) analysed.', numel(app.kine)));
        onFishSelected();
    end

    % ---- Animal selected ----
    function onFishSelected()
        if isempty(app.kine), return; end
        names = {app.kine.name};
        fi    = find(strcmp(app.fishList.Value, names), 1);
        if isempty(fi), return; end

        k  = app.kine(fi);
        fp = app.fp(fi);
        has_z = isfield(fp, 'has_z') && fp.has_z;

        % -- Results text --
        L = {};
        L{end+1} = sprintf('Animal:             %s', k.name);
        L{end+1} = repmat('-',1,34);
        L{end+1} = sprintf('Head TBF:           %.3f Hz',  k.head_TBF);
        L{end+1} = sprintf('Tail TBF:           %.3f Hz',  k.tail_TBF);
        if has_z
            L{end+1} = sprintf('Head TBF (Z):       %.3f Hz',  k.headZ_TBF);
            L{end+1} = sprintf('Tail TBF (Z):       %.3f Hz',  k.tailZ_TBF);
        end
        L{end+1} = repmat('-',1,34);
        L{end+1} = sprintf('Head amp (Y):       %.4f BL',  k.headAmp);
        L{end+1} = sprintf('Tail amp (Y):       %.4f BL',  k.tailAmp);
        L{end+1} = sprintf('Head/tail ratio:    %.4f',     k.headTailAmpRatio);
        L{end+1} = sprintf('Min amp (Y):        %.4f BL @ %.3f', k.minAmp, k.minAmpLoc);
        L{end+1} = sprintf('Max amp (Y):        %.4f BL @ %.3f', k.maxAmp, k.maxAmpLoc);
        if has_z
            L{end+1} = repmat('-',1,34);
            L{end+1} = sprintf('Head amp (Z):       %.4f',     k.headAmpZ);
            L{end+1} = sprintf('Tail amp (Z):       %.4f',     k.tailAmpZ);
            L{end+1} = sprintf('Min amp (Z):        %.4f @ %.3f', k.minAmpZ, k.minAmpZLoc);
            L{end+1} = sprintf('Max amp (Z):        %.4f @ %.3f', k.maxAmpZ, k.maxAmpZLoc);
        end
        L{end+1} = repmat('-',1,34);
        L{end+1} = sprintf('Wavelength:         %.4f BL',  k.wavelength);
        L{end+1} = repmat('-',1,34);
        L{end+1} = sprintf('Max curv (XY):      %.4f @ %.3f', k.maxCurv, k.maxCurvLoc);
        if has_z
            L{end+1} = sprintf('Max curv (3D):      %.4f @ %.3f', k.maxCurv3D, k.maxCurv3DLoc);
        end
        app.resultsArea.Value = L;

        % -- Midlines plot --
        drawMidlines(k, fp);

        % -- Plot 2: Amplitude envelope --
        cla(app.ax(2)); hold(app.ax(2),'on');
        s = k.s_norm;
        fill(app.ax(2), [s fliplr(s)], [k.amp_mean+k.amp_std fliplr(k.amp_mean-k.amp_std)], ...
             [0.6 0.8 1], 'EdgeColor','none','FaceAlpha',0.4);
        plot(app.ax(2), s, k.amp_mean, 'b-', 'LineWidth', 2, 'DisplayName','Y amp');
        if has_z && ~isempty(k.ampZ_mean)
            fill(app.ax(2), [s fliplr(s)], [k.ampZ_mean+k.ampZ_std fliplr(k.ampZ_mean-k.ampZ_std)], ...
                 [1 0.75 0.6], 'EdgeColor','none','FaceAlpha',0.3);
            plot(app.ax(2), s, k.ampZ_mean, 'r-', 'LineWidth', 2, 'DisplayName','Z amp');
        end
        xline(app.ax(2), k.minAmpLoc, 'b:'); xline(app.ax(2), k.maxAmpLoc, 'b--');
        legend(app.ax(2), 'Location','northwest','FontSize',9);
        hold(app.ax(2),'off'); xlim(app.ax(2),[0 1]);

        % -- Plot 3: Curvature --
        cla(app.ax(3)); hold(app.ax(3),'on');
        fill(app.ax(3), [s fliplr(s)], [k.curv_mean+k.curv_std fliplr(k.curv_mean-k.curv_std)], ...
             [1 0.8 0.7], 'EdgeColor','none','FaceAlpha',0.4);
        plot(app.ax(3), s, k.curv_mean, 'r-', 'LineWidth', 2, 'DisplayName','XY curv');
        if has_z && ~isempty(k.curv3d_mean)
            plot(app.ax(3), s, k.curv3d_mean, 'm--', 'LineWidth', 1.5, 'DisplayName','3D curv');
        end
        xline(app.ax(3), k.maxCurvLoc, 'k--');
        legend(app.ax(3), 'Location','northwest','FontSize',9);
        hold(app.ax(3),'off'); xlim(app.ax(3),[0 1]);

        % -- Plot 4: FFT spectra --
        % NOTE: Unclean FFT spectra are expected — they reflect real broadband
        % signal content in the time series (e.g. non-sinusoidal waveforms,
        % aperiodic motion, NaN-gap interpolation artefacts). This is a data
        % characteristic, not a coding bug. The dominant_freq() function uses
        % min_freq masking to reliably find the peak regardless.
        cla(app.ax(4)); hold(app.ax(4),'on');
        plot(app.ax(4), k.head_fft_freq, k.head_fft_power, 'b-', 'LineWidth',1.2, ...
             'DisplayName', sprintf('Head Y (%.2fHz)', k.head_TBF));
        plot(app.ax(4), k.tail_fft_freq, k.tail_fft_power, 'r-', 'LineWidth',1.2, ...
             'DisplayName', sprintf('Tail Y (%.2fHz)', k.tail_TBF));
        minFreq = app.minFreqField.Value;
        xline(app.ax(4), minFreq,     'k:', 'LineWidth',1, ...
              'Label', sprintf('min %.1fHz',minFreq));
        xline(app.ax(4), k.head_TBF, 'b--', 'LineWidth',1);
        xline(app.ax(4), k.tail_TBF, 'r--', 'LineWidth',1);
        legend(app.ax(4),'Location','northeast','FontSize',9);
        hold(app.ax(4),'off');
    end

    % ---- Draw midlines (called on fish select AND on control change) ----
    function drawMidlines(k, fp)
        if nargin < 2
            % Re-draw using current selection
            if isempty(app.kine), return; end
            names = {app.kine.name};
            fi = find(strcmp(app.fishList.Value, names), 1);
            if isempty(fi), return; end
            k  = app.kine(fi);
            fp = app.fp(fi);
        end

        cla(app.ax(1));
        valid = find(~any(isnan(k.Y_interp), 2));
        if isempty(valid), return; end

        nShow  = min(app.midlineCount.Value, numel(valid));
        % Evenly spaced indices across valid frames
        idx_show = round(linspace(1, numel(valid), nShow));
        idx_show = unique(idx_show);

        % Thickness: thicker when fewer lines shown (range 0.5–3.5)
        lw = max(0.5, min(3.5, 0.5 + 30 / max(nShow, 1)));

        % Alpha: slightly more opaque when fewer
        alpha_val = max(0.15, min(0.85, 0.15 + 10 / max(nShow, 1)));

        cmapName = app.cmapDrop.Value;
        cmap = feval(cmapName, numel(idx_show));

        hold(app.ax(1),'on');
        yline(app.ax(1), 0, 'k--', 'LineWidth', 0.8);

        for ii = 1:numel(idx_show)
            f = valid(idx_show(ii));
            plot(app.ax(1), k.X_interp(f,:), k.Y_interp(f,:), ...
                 'Color', [cmap(ii,:) alpha_val], 'LineWidth', lw);
        end

        % Raw measured points at middle valid frame
        mf = valid(round(end/2));
        scatter(app.ax(1), fp.X(mf,:), fp.Y(mf,:), 50, 'r', 'filled', ...
                'DisplayName','Measured pts (mid frame)', 'HandleVisibility','on');

        % Colorbar showing frame number
        colormap(app.ax(1), cmap);
        cb = colorbar(app.ax(1));
        cb.Label.String = 'Frame (early -> late)';
        clim(app.ax(1), [valid(idx_show(1)), valid(idx_show(end))]);

        % Y limits: use actual data range + 20% padding (NOT axis equal —
        % the fish body is much longer than it is wide, so equal squashes the motion)
        all_y = k.Y_interp(~isnan(k.Y_interp));
        if ~isempty(all_y)
            y_pad = (max(all_y) - min(all_y)) * 0.20 + 0.01;
            ylim(app.ax(1), [min(all_y) - y_pad,  max(all_y) + y_pad]);
        end
        % X: tight to data
        all_x = k.X_interp(~isnan(k.X_interp));
        if ~isempty(all_x)
            xlim(app.ax(1), [min(all_x), max(all_x)]);
        end

        hold(app.ax(1),'off');
        % Count how many frames were skipped (NaN in interpolated midline)
        nSkipped = sum(any(isnan(k.Y_interp), 2));
        title(app.ax(1), sprintf('Midlines (n=%d shown, %d frames skipped — occluded middle points)', ...
              numel(idx_show), nSkipped), 'FontSize', 9);
    end

    function refreshMidlines()
        if isempty(app.kine), return; end
        drawMidlines();
    end

    %% ================================================================
    %  FIN ANALYSIS  (computation delegated to compute_fin_kinematics.m)
    %% ================================================================

    function onComputeFin()
        if isempty(app.fp)
            uialert(fig, 'Run kinematics analysis first.', 'No data'); return
        end
        if ~isfield(app.fp, 'has_z') || ~app.fp(1).has_z
            uialert(fig, '3D data required for fin analysis. Load a file with X, Y, Z columns.', '3D required');
            return
        end

        rootName = app.finRootDrop.Value;
        tipName  = app.finTipDrop.Value;

        if strcmp(rootName, tipName)
            uialert(fig, 'Root and tip must be different points.', 'Selection'); return
        end
        if contains(rootName, 'N/A') || contains(tipName, 'N/A')
            uialert(fig, 'Named-format file required for fin analysis.', 'Format'); return
        end

        csvPath = strtrim(app.fileField.Value);
        fps     = app.fpsField.Value;

        setStatus('Computing fin kinematics...');
        try
            fd = compute_fin_kinematics(csvPath, rootName, tipName, fps);
        catch ME
            uialert(fig, ME.message, 'Fin computation error');
            setStatus('Fin computation failed.'); return
        end
        app.fin_data = fd;
        setStatus('Fin kinematics done.');

        valid  = fd.valid;
        frames = fd.frames;

        % ---- Results text ----
        R = {};
        R{end+1} = sprintf('Fin:  root = %s   |   tip = %s', rootName, tipName);
        R{end+1} = repmat('-',1,60);
        R{end+1} = sprintf('Fin length:          %.4f +/- %.4f (mean +/- SD)', fd.mean_length, fd.std_length);
        R{end+1} = sprintf('Total tip distance:  %.4f', fd.total_dist);
        R{end+1} = sprintf('Mean tip speed:      %.4f /s  (SD %.4f  peak %.4f)', ...
                            fd.mean_speed, fd.std_speed, fd.peak_speed);
        R{end+1} = repmat('-',1,60);
        R{end+1} = sprintf('Yaw   — mean: %7.2f deg   SD: %6.2f   range: %.2f', ...
                            fd.mean_yaw,   fd.std_yaw,   fd.range_yaw);
        R{end+1} = sprintf('Pitch — mean: %7.2f deg   SD: %6.2f   range: %.2f', ...
                            fd.mean_pitch, fd.std_pitch, fd.range_pitch);
        R{end+1} = sprintf('Roll  — mean: %7.2f deg   SD: %6.2f   range: %.2f', ...
                            fd.mean_roll,  fd.std_roll,  fd.range_roll);
        R{end+1} = '  (Roll = atan2(vz,vy): fin cant in body cross-section plane)';
        R{end+1} = repmat('-',1,60);
        R{end+1} = sprintf('Mean angular velocity: %.4f deg/s  (SD %.4f  peak %.4f)', ...
                            fd.mean_ang_vel, fd.std_ang_vel, fd.peak_ang_vel);
        R{end+1} = sprintf('Valid frames:          %d / %d (%.1f%%)', ...
                            fd.n_valid, fd.n_frames, fd.pct_valid);
        app.finResultsArea.Value = R;

        % ---- Plot 1: Pitch / Roll / Yaw over time ----
        cla(app.finAx(1)); hold(app.finAx(1),'on');
        plot(app.finAx(1), frames(valid), fd.pitch(valid), 'b-',  'LineWidth',1.4, 'DisplayName','Pitch');
        plot(app.finAx(1), frames(valid), fd.yaw(valid),   'r-',  'LineWidth',1.4, 'DisplayName','Yaw');
        plot(app.finAx(1), frames(valid), fd.roll(valid),  'Color',[0.1 0.7 0.3], ...
             'LineWidth',1.4, 'DisplayName','Roll');
        yline(app.finAx(1), 0, 'k--', 'LineWidth',0.8, 'HandleVisibility','off');
        legend(app.finAx(1), 'Location','best','FontSize',9);
        xlabel(app.finAx(1),'Frame'); ylabel(app.finAx(1),'Angle (deg)');
        title(app.finAx(1),'Pitch / Yaw / Roll over time','FontSize',10);
        hold(app.finAx(1),'off');

        % ---- Plot 2: Cumulative distance ----
        cla(app.finAx(2));
        plot(app.finAx(2), frames(valid), fd.cum_dist(valid), 'Color',[0.1 0.55 0.3], 'LineWidth',1.8);
        xlabel(app.finAx(2),'Frame'); ylabel(app.finAx(2),'Cumulative distance');
        title(app.finAx(2), sprintf('Tip distance (total = %.3f)', fd.total_dist), 'FontSize',10);

        % ---- Plot 3: Tip trajectory in XY — connected line colored by frame ----
        cla(app.finAx(3));
        tv    = fd.tip_xyz(valid,:);
        fv    = frames(valid);
        nv    = size(tv,1);
        cmap_xy = parula(max(nv,2));
        hold(app.finAx(3),'on');
        for seg = 2:nv
            plot(app.finAx(3), tv(seg-1:seg,1), tv(seg-1:seg,2), '-', ...
                 'Color', cmap_xy(seg,:), 'LineWidth', 2.0, 'HandleVisibility','off');
        end
        scatter(app.finAx(3), tv(1,1),   tv(1,2),   60, 'g', 'filled', ...
                'MarkerEdgeColor','k', 'DisplayName','Start');
        scatter(app.finAx(3), tv(end,1), tv(end,2), 60, 'r', 'filled', ...
                'MarkerEdgeColor','k', 'DisplayName','End');
        colormap(app.finAx(3), parula);
        cb_xy = colorbar(app.finAx(3));
        cb_xy.Label.String = 'Frame';
        clim(app.finAx(3), [fv(1) fv(end)]);
        legend(app.finAx(3),'Location','best','FontSize',8);
        xlabel(app.finAx(3),'X'); ylabel(app.finAx(3),'Y');
        title(app.finAx(3),'Tip XY trajectory (colored by frame)','FontSize',10);
        hold(app.finAx(3),'off');

        % ---- Plot 4: 3D tip trajectory ----
        cla(app.finAx3D);
        nValid = size(tv,1);
        cmap3d = parula(max(nValid, 2));
        hold(app.finAx3D,'on');
        for seg = 2:nValid
            plot3(app.finAx3D, tv(seg-1:seg,1), tv(seg-1:seg,2), tv(seg-1:seg,3), ...
                  '-', 'Color',[cmap3d(seg,:) 0.85], 'LineWidth',2.5, ...
                  'HandleVisibility','off');
        end
        scatter3(app.finAx3D, tv(1,1),   tv(1,2),   tv(1,3),   90, 'g', 'filled', ...
                 'DisplayName','Start', 'MarkerEdgeColor','w');
        scatter3(app.finAx3D, tv(end,1), tv(end,2), tv(end,3), 90, 'r', 'filled', ...
                 'DisplayName','End',   'MarkerEdgeColor','w');
        legend(app.finAx3D,'Location','best','FontSize',9,'TextColor',[1 1 1],'Color',[0.1 0.1 0.2]);
        xlabel(app.finAx3D,'X','Color',[0.8 0.8 0.8],'FontSize',10);
        ylabel(app.finAx3D,'Y','Color',[0.8 0.8 0.8],'FontSize',10);
        zlabel(app.finAx3D,'Z','Color',[0.8 0.8 0.8],'FontSize',10);
        title(app.finAx3D, sprintf('Fin tip trajectory (3D)\n%s  \x2192  %s', rootName, tipName), ...
              'FontSize', 10, 'Color', [0.95 0.95 0.5], 'FontWeight','bold');
        colormap(app.finAx3D, parula);
        cb3 = colorbar(app.finAx3D);
        cb3.Label.String = 'Frame';
        cb3.Color = [0.8 0.8 0.8];
        cb3.Label.Color = [0.8 0.8 0.8];
        hold(app.finAx3D,'off');
        view(app.finAx3D, 35, 25);
        axis(app.finAx3D,'tight');
        xl=xlim(app.finAx3D); yl=ylim(app.finAx3D); zl=zlim(app.finAx3D);
        pad3 = max([diff(xl) diff(yl) diff(zl)]) * 0.12 + 0.01;
        xlim(app.finAx3D, xl+[-pad3 pad3]);
        ylim(app.finAx3D, yl+[-pad3 pad3]);
        zlim(app.finAx3D, zl+[-pad3 pad3]);
    end

    % ---- Fin animation in a new figure ----
    function onAnimateFin()
        if isempty(app.fin_data)
            uialert(fig, 'Compute fin analysis first.', 'No data'); return
        end
        runFinAnimation(app.fin_data);
    end

    function runFinAnimation(fd)
        valid    = fd.valid;
        frames_v = find(valid);
        if numel(frames_v) < 2
            uialert(fig,'Not enough valid frames to animate.','Animation'); return
        end

        vx_v = fd.vx(frames_v);
        vy_v = fd.vy(frames_v);
        vz_v = fd.vz(frames_v);

        animFig = figure('Name', sprintf('Fin Animation: %s -> %s', fd.rootName, fd.tipName), ...
                         'Position', [150 80 900 740], 'Color', [0.05 0.05 0.1]);

        % ---- Control bar at bottom ----
        uicontrol(animFig, 'Style','pushbutton', 'String','▶  Restart', ...
            'Units','pixels', 'Position',[20 10 120 30], ...
            'BackgroundColor',[0.18 0.42 0.75], 'ForegroundColor',[1 1 1], ...
            'FontSize',11, 'FontWeight','bold', ...
            'Callback', @(~,~) runFinAnimation(fd));   % re-call this function

        uicontrol(animFig, 'Style','pushbutton', 'String','✕  Close', ...
            'Units','pixels', 'Position',[155 10 100 30], ...
            'BackgroundColor',[0.6 0.15 0.15], 'ForegroundColor',[1 1 1], ...
            'FontSize',11, ...
            'Callback', @(~,~) close(animFig));

        speed_lbl = uicontrol(animFig, 'Style','text', 'String','Speed:', ...
            'Units','pixels', 'Position',[275 12 48 22], ...
            'BackgroundColor',[0.05 0.05 0.1], 'ForegroundColor',[0.85 0.85 0.85], ...
            'FontSize',10, 'HorizontalAlignment','right'); %#ok<NASGU>
        speed_slider = uicontrol(animFig, 'Style','slider', ...
            'Units','pixels', 'Position',[328 14 160 18], ...
            'Min',0.1, 'Max',8, 'Value',1, ...
            'SliderStep',[0.05 0.2]);
        uicontrol(animFig, 'Style','text', 'String','0.1x', ...
            'Units','pixels', 'Position',[275 11 50 18], ...
            'BackgroundColor',[0.05 0.05 0.1], 'ForegroundColor',[0.6 0.6 0.6],'FontSize',9);
        uicontrol(animFig, 'Style','text', 'String','8x', ...
            'Units','pixels', 'Position',[490 11 30 18], ...
            'BackgroundColor',[0.05 0.05 0.1], 'ForegroundColor',[0.6 0.6 0.6],'FontSize',9);

        % Axes — leave room for control bar
        ax_anim = axes(animFig, ...
            'Units','pixels', 'Position',[60 55 800 660], ...
            'Color',  [0.05 0.05 0.1], ...
            'XColor', [0.75 0.75 0.75], ...
            'YColor', [0.75 0.75 0.75], ...
            'ZColor', [0.75 0.75 0.75], ...
            'GridColor', [0.3 0.3 0.3], 'GridAlpha', 0.4, 'FontSize', 10);
        hold(ax_anim,'on'); grid(ax_anim,'on');
        xlabel(ax_anim,'X','FontSize',12,'Color',[0.85 0.85 0.85]);
        ylabel(ax_anim,'Y','FontSize',12,'Color',[0.85 0.85 0.85]);
        zlabel(ax_anim,'Z','FontSize',12,'Color',[0.85 0.85 0.85]);
        title(ax_anim, sprintf('Fin vector:  %s  \x2192  %s\n(root fixed at origin)', ...
              fd.rootName, fd.tipName), ...
              'Color',[0.95 0.95 0.5],'FontSize',12,'FontWeight','bold');
        view(ax_anim, 35, 25);

        % Axis limits from vector data
        vrange = max([range(vx_v) range(vy_v) range(vz_v)]);
        pad_a  = max(vrange * 0.28, 0.05);
        xlim(ax_anim, [min(vx_v)-pad_a  max(vx_v)+pad_a]);
        ylim(ax_anim, [min(vy_v)-pad_a  max(vy_v)+pad_a]);
        zlim(ax_anim, [min(vz_v)-pad_a  max(vz_v)+pad_a]);

        % Ghost trail
        plot3(ax_anim, vx_v, vy_v, vz_v, '-', ...
              'Color',[0.6 0.6 0.6 0.2],'LineWidth',0.8,'HandleVisibility','off');

        % Growing trail
        h_trail = plot3(ax_anim, vx_v(1), vy_v(1), vz_v(1), '-', ...
                        'Color',[0.4 0.75 1.0 0.6],'LineWidth',2.0,'HandleVisibility','off');

        % Root
        scatter3(ax_anim, 0,0,0, 130,'g','filled','MarkerEdgeColor','w','DisplayName','Root (origin)');

        % Fin vector line
        h_vec = plot3(ax_anim,[0 vx_v(1)],[0 vy_v(1)],[0 vz_v(1)], ...
                      'w-','LineWidth',3.5,'HandleVisibility','off');

        % Tip sphere
        h_tip = scatter3(ax_anim, vx_v(1),vy_v(1),vz_v(1), ...
                         240,'y','filled','MarkerEdgeColor','w','DisplayName','Fin tip');

        legend(ax_anim,'Location','northeast','TextColor','white','Color',[0.1 0.1 0.2],'FontSize',10);

        % HUD overlays
        h_frame = text(ax_anim,0.02,0.97,0.02, sprintf('Frame %d / %d',frames_v(1),frames_v(end)), ...
                       'Units','normalized','Color',[1 1 1],'FontSize',12,'FontWeight','bold', ...
                       'VerticalAlignment','top');
        h_ang   = text(ax_anim,0.02,0.90,0.02,'', ...
                       'Units','normalized','Color',[1 0.95 0.4],'FontSize',11, ...
                       'VerticalAlignment','top');

        % Animation loop
        target_fps = 25;
        nF = numel(frames_v);
        for ii = 1:nF
            if ~isvalid(animFig), return; end
            f = frames_v(ii);

            h_vec.XData = [0 vx_v(ii)];
            h_vec.YData = [0 vy_v(ii)];
            h_vec.ZData = [0 vz_v(ii)];
            h_tip.XData = vx_v(ii);
            h_tip.YData = vy_v(ii);
            h_tip.ZData = vz_v(ii);
            h_trail.XData = vx_v(1:ii);
            h_trail.YData = vy_v(1:ii);
            h_trail.ZData = vz_v(1:ii);

            h_frame.String = sprintf('Frame %d / %d', f, frames_v(end));
            h_ang.String   = sprintf('Yaw: %.1f\x00B0   Pitch: %.1f\x00B0   Roll: %.1f\x00B0', ...
                                     fd.yaw(f), fd.pitch(f), fd.roll(f));

            drawnow limitrate;
            spd = speed_slider.Value;
            pause(1 / (target_fps * spd));
        end
    end

    %% ================================================================
    %  HELPERS
    %% ================================================================

    function setStatus(msg)
        app.statusLabel.Text = msg; drawnow;
    end

    function y_out = section_label(parent, txt, y_in)
        uilabel(parent, 'Text', txt, 'Position', [12 y_in-18 316 18], ...
                'FontSize', 10, 'FontWeight', 'bold', 'FontColor', [0.5 0.5 0.5]);
        y_out = y_in - 24;
    end

end