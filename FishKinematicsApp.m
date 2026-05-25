function FishKinematicsApp()
% FISHKINEMATICSAPP  Interactive UI for fish midline kinematics + fin analysis.
%
%   Supports four CSV formats:
%     FORMAT A (DLC-style):      columns Fish1_P1_x, Fish1_P1_y[, Fish1_P1_z] ...
%     FORMAT B (named):          columns eye_X, eye_Y[, eye_Z], snout_X ...
%     FORMAT C (numbered 3D/2D): columns pt1_X, pt1_Y[, pt1_Z], pt2_X, ... (12 anatomical pts)
%     FORMAT D (dual-camera 2D): columns pt1_cam1_X, pt1_cam1_Y, pt1_cam2_X, ... (cam1=lateral, cam2=ventral)
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
    app.pect_data = [];   % pectoral phase result (Format D)
    % Axis mapping: which signed CSV axis maps to each world role.
    % Default = identity (CSV X→top-X, CSV Y→top-Y, CSV Z→vertical).
    app.axis_map.top_x = '+X';
    app.axis_map.top_y = '+Y';
    app.axis_map.vert  = '+Z';
    % CSV export state
    app.csv_rows = {};    % cell array of flat structs, one per exported analysis
    app.csv_path = '';    % last used output file path

    %% ---- Figure ----
    fig = uifigure('Name',     'Kinemetrix', ...
                   'Position', [60 40 1260 860], ...
                   'Color',    [0.95 0.95 0.95]);

    % ---- Title bar ----
    titleBar = uipanel(fig, 'Position', [0 820 1260 40], ...
        'BackgroundColor', [0.12 0.22 0.45], 'BorderType', 'none');
    uilabel(titleBar, 'Text', 'Kinemetrix', ...
        'Position', [0 0 1260 40], ...
        'FontSize', 22, 'FontWeight', 'bold', ...
        'FontColor', [1 1 1], ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment',   'center');

    %% ================================================================
    %  LEFT PANEL  (inputs, point selector, results)
    %% ================================================================
    LP = uipanel(fig, 'Position', [10 10 340 806], ...
                 'BackgroundColor', [1 1 1], 'BorderType', 'line', ...
                 'Title', '', 'Scrollable', 'on');

    y = 1100;   % tall virtual canvas — LP scrolls to show all content

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

    % ---- Axis Mapping (shown only for 3D formats) ----
    axMapH = 122;
    app.axMapPanel = uipanel(LP, 'Position', [8 y-axMapH 324 axMapH], ...
        'BackgroundColor', [1 0.97 0.94], 'BorderType', 'line', 'Title', 'Axis Mapping (3D)');
    app.axMapPanel.Visible = 'off';

    axis_opts = {'+X','-X','+Y','-Y','+Z','-Z'};

    uilabel(app.axMapPanel, 'Text', 'Top-view  X  →  CSV axis:', ...
        'Position', [6 76 168 18], 'FontSize', 10);
    app.axTopXDrop = uidropdown(app.axMapPanel, ...
        'Items', axis_opts, 'Value', '+X', ...
        'Position', [178 74 126 22], 'FontSize', 10);

    uilabel(app.axMapPanel, 'Text', 'Top-view  Y  →  CSV axis:', ...
        'Position', [6 50 168 18], 'FontSize', 10);
    app.axTopYDrop = uidropdown(app.axMapPanel, ...
        'Items', axis_opts, 'Value', '+Y', ...
        'Position', [178 48 126 22], 'FontSize', 10);

    uilabel(app.axMapPanel, 'Text', 'Vertical   Z  →  CSV axis:', ...
        'Position', [6 24 168 18], 'FontSize', 10);
    app.axVertDrop = uidropdown(app.axMapPanel, ...
        'Items', axis_opts, 'Value', '+Z', ...
        'Position', [178 22 126 22], 'FontSize', 10);

    app.axWarnLabel = uilabel(app.axMapPanel, ...
        'Text', '', 'Position', [6 3 308 16], ...
        'FontSize', 9, 'FontColor', [0.75 0.2 0.1], 'WordWrap', 'off');

    app.axTopXDrop.ValueChangedFcn = @(~,~) validateAxisMap();
    app.axTopYDrop.ValueChangedFcn = @(~,~) validateAxisMap();
    app.axVertDrop.ValueChangedFcn  = @(~,~) validateAxisMap();

    y = y - (axMapH + 8);

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

    uilabel(app.ptPanel, 'Text', 'Selected order:', ...
        'Position', [216 162 110 18], 'FontSize', 10, 'FontWeight', 'bold', 'FontColor', [0.4 0.4 0.4]);
    app.selList = uilistbox(app.ptPanel, 'Position', [216 20 100 142], ...
        'Multiselect', 'on', 'Items', {});
    app.selHintLabel = uilabel(app.ptPanel, ...
        'Text', '↑ head  /  tail ↓', ...
        'Position', [216 3 100 16], 'FontSize', 9, ...
        'HorizontalAlignment', 'center', 'FontColor', [0.45 0.45 0.65]);

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
    y = y - 28;

    % ---- Export CSV ----
    y = section_label(LP, 'EXPORT', y);
    app.csvRowLabel = uilabel(LP, 'Text', '0 rows staged', ...
        'Position', [12 y-18 180 18], 'FontSize', 10, 'FontColor', [0.35 0.35 0.35]);
    uibutton(LP, 'Text', 'Clear staged', ...
        'Position', [196 y-20 132 20], 'FontSize', 9, ...
        'BackgroundColor', [0.92 0.88 0.88], ...
        'ButtonPushedFcn', @(~,~) onClearCSV());
    y = y - 26;
    uibutton(LP, 'Text', '⬇  Export kinematics CSV', ...
        'Position', [12 y-30 316 30], 'FontSize', 11, 'FontWeight', 'bold', ...
        'BackgroundColor', [0.13 0.52 0.30], 'FontColor', [1 1 1], ...
        'ButtonPushedFcn', @(~,~) onExportCSV());
    y = y - 38;

    % ---- Animal selector ----
    y = section_label(LP, 'SELECT ANIMAL', y);
    app.fishList = uilistbox(LP, 'Position', [12 y-70 316 68], ...
        'Items', {'(run analysis first)'}, 'ValueChangedFcn', @(~,~) onFishSelected());
    y = y - 78;

    % ---- Results ----
    y = section_label(LP, 'KINEMATIC VALUES', y);
    app.resultsArea = uitextarea(LP, 'Position', [12 y-220 316 214], ...
        'Editable', 'off', 'FontSize', 11, 'FontName', 'Courier New', ...
        'BackgroundColor', [0.97 0.97 0.97], 'Value', {''});
    y = y - 228;

    %% ================================================================
    %  RIGHT PANEL  — Tabbed: Kinematics | Fin Analysis
    %% ================================================================
    tg = uitabgroup(fig, 'Position', [360 10 880 806]);

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
    xlabel(app.ax(1),'Body position (BL,  0=head → 1=tail)');   ylabel(app.ax(1),'Lateral displacement (BL)');
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

    % ---- Tab 3: Pectoral Fin Phase (Format D dual-camera) ----
    tabPect = uitab(tg, 'Title', 'Pectoral Phase (2D)');

    app.pectInfoLabel = uilabel(tabPect, 'Position', [8 756 860 36], ...
        'Text', 'Load a dual-camera CSV (Format D: pt1_cam1_X ...) and run analysis to see pectoral fin phase.', ...
        'FontSize', 12, 'WordWrap', 'on', 'FontColor', [0.4 0.4 0.6]);

    app.pectResultArea = uitextarea(tabPect, 'Position', [8 580 860 168], ...
        'Editable', 'off', 'FontSize', 12, 'FontName', 'Courier New', ...
        'BackgroundColor', [0.97 0.97 0.97]);

    % Two axes: time-domain overlay and cross-correlation
    app.pectAx = gobjects(2,1);
    app.pectAx(1) = uiaxes(tabPect, 'Position', [8 300 860 272], ...
        'BackgroundColor', [1 1 1], 'Box', 'on');
    title(app.pectAx(1), 'Pectoral fin tip signals over time (right pt2 vs left pt12)', 'FontSize', 10);
    xlabel(app.pectAx(1), 'Frame');  ylabel(app.pectAx(1), 'Y position');
    app.pectAx(1).Toolbar.Visible = 'off';
    text(app.pectAx(1), 0.5, 0.5, 'No data', 'Units','normalized', ...
         'HorizontalAlignment','center','FontSize',12,'Color',[0.75 0.75 0.75]);

    app.pectAx(2) = uiaxes(tabPect, 'Position', [8 8 860 280], ...
        'BackgroundColor', [1 1 1], 'Box', 'on');
    title(app.pectAx(2), 'Cross-correlation (right pectoral vs left pectoral)', 'FontSize', 10);
    xlabel(app.pectAx(2), 'Lag (frames)');  ylabel(app.pectAx(2), 'Normalized cross-correlation');
    app.pectAx(2).Toolbar.Visible = 'off';
    text(app.pectAx(2), 0.5, 0.5, 'No data', 'Units','normalized', ...
         'HorizontalAlignment','center','FontSize',12,'Color',[0.75 0.75 0.75]);

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
        % Allow multi-select so CURVES files (one per fish) can be loaded together
        [fn, fp_] = uigetfile('*.csv', 'Select tracking CSV (hold Ctrl/Shift for multiple CURVES files)', ...
                              'MultiSelect', 'on');
        if isequal(fn, 0), return; end
        if ischar(fn)
            % Single file
            fullpath = fullfile(fp_, fn);
            app.fileField.Value = fullpath;
            detectFormat(fullpath);
        else
            % Multiple files — assume CURVES format
            app.fmt = 'curves';
            % Store all paths as semicolon-separated in fileField
            paths = cellfun(@(f) fullfile(fp_, f), fn, 'UniformOutput', false);
            app.fileField.Value = strjoin(paths, ';');
            n_fish = numel(fn);
            app.fmtLabel.Text = sprintf('Format E: CURVES midline  (%d files selected)', n_fish);
            app.ptPanel.Visible = 'off';
            app.axMapPanel.Visible = 'off';
            app.finRootDrop.Items = {'N/A — 2D CURVES format'};
            app.finTipDrop.Items  = {'N/A — 2D CURVES format'};
        end
    end

    function detectFormat(path)
        opts = detectImportOptions(path);
        opts.VariableNamingRule = 'preserve';
        opts.DataLines = [1 2];
        T_peek = readtable(path, opts);
        cols = T_peek.Properties.VariableNames;

        % Helper: show or hide the axis mapping panel
        function showAxisPanel(tf)
            if tf
                app.axMapPanel.Visible = 'on';
                validateAxisMap();
            else
                app.axMapPanel.Visible = 'off';
            end
        end

        % FORMAT E: CURVES — first col name is a number (body station in mm), second is blank/auto-named.
        % With VariableNamingRule='preserve', MATLAB names empty CSV columns as 'Var2','Var4',... 
        % (not 'Unnamed...'), so we must accept both patterns.
        is_col2_filler = numel(cols) > 1 && ( ...
            isempty(cols{2}) || ...
            ~isempty(regexp(cols{2}, '^Unnamed', 'once')) || ...
            ~isempty(regexp(cols{2}, '^Var\d+$',  'once')) );
        is_curves = ~isempty(regexp(cols{1}, '^\d+$', 'once')) && is_col2_filler;
        if is_curves
            num_cols = cols(~cellfun(@isempty, regexp(cols, '^\d+$')));
            app.fmt = 'curves';
            app.fmtLabel.Text = sprintf('Format E: CURVES midline  (%d body stations, pre-transformed X/Y in BL)', numel(num_cols));
            app.ptPanel.Visible = 'off';
            showAxisPanel(false);   % 2D pre-transformed — axis mapping not applicable
            app.finRootDrop.Items = {'N/A — 2D CURVES format'};
            app.finTipDrop.Items  = {'N/A — 2D CURVES format'};
            return;
        end

        % FORMAT D: pt<N>_cam<M>_X — dual-camera 2D
        if any(~cellfun(@isempty, regexp(cols, '^pt\d+_cam\d+_[XxYy]$')))
            app.fmt = 'dual_camera';
            app.fmtLabel.Text = 'Format D: Dual-camera 2D  (pt1_cam1_X ... | cam1=lateral, cam2=ventral)';
            app.ptPanel.Visible = 'off';
            showAxisPanel(false);   % 2D only
            app.finRootDrop.Items = {'N/A — 2D dual-camera format'};
            app.finTipDrop.Items  = {'N/A — 2D dual-camera format'};
            return;
        end

        % FORMAT C: pt<N>_X — numbered 3D or 2D points
        if any(~cellfun(@isempty, regexp(cols, '^pt\d+_[XxYyZz]$')))
            app.fmt = 'numbered';
            has_z_c = any(~cellfun(@isempty, regexp(cols, '^pt\d+_[Zz]$')));
            dimStr  = '3D';
            if ~has_z_c, dimStr = '2D'; end
            app.fmtLabel.Text = sprintf('Format C: Numbered points %s  (pt1_X, pt2_X ... | 12 anatomical pts)', dimStr);
            app.ptPanel.Visible = 'off';
            showAxisPanel(has_z_c);
            % Populate fin dropdowns from numbered point anatomy
            tok_x   = regexp(cols, '^(pt\d+)_[Xx]$', 'tokens');
            pt_b    = cellfun(@(t) t{1}{1}, tok_x(~cellfun(@isempty,tok_x)), 'UniformOutput', false);
            pt_b    = unique(pt_b, 'stable');
            pt_nums = cellfun(@(s) str2double(regexp(s,'\d+','match','once')), pt_b);
            [~, ord] = sort(pt_nums);
            pt_b    = pt_b(ord);
            app.finRootDrop.Items = pt_b;
            app.finTipDrop.Items  = pt_b;
            % Suggest anatomical defaults: pt1=pect_base_R, pt2=pect_tip_R
            if numel(pt_b) >= 2
                app.finRootDrop.Value = pt_b{1};   % pt1 = pect_base_R
                app.finTipDrop.Value  = pt_b{2};   % pt2 = pect_tip_R
            end
            return;
        end

        % FORMAT A: DLC FishN_PN_x
        dlc_match = ~cellfun(@isempty, regexp(cols, '^Fish\d+_P\d+_[xyXY]'));
        if any(dlc_match)
            has_z_dlc = any(~cellfun(@isempty, regexp(cols, '^Fish\d+_P\d+_[zZ]$')));
            app.fmt = 'DLC';
            app.fmtLabel.Text = 'Format A: DLC  (Fish1_P1_x columns)';
            app.ptPanel.Visible = 'off';
            showAxisPanel(has_z_dlc);
            app.finRootDrop.Items = {'N/A — use named format'};
            app.finTipDrop.Items  = {'N/A — use named format'};
        else
            % FORMAT B: named columns
            has_z_b = any(~cellfun(@isempty, regexp(cols, '^.+_[Zz]$')));
            app.fmt = 'named';
            app.fmtLabel.Text = 'Format B: Named  (eye_X, snout_Y ... columns)';
            app.ptPanel.Visible = 'on';
            showAxisPanel(has_z_b);
            tok = regexp(cols, '^(.+)_[Xx]$', 'tokens');
            bases = cellfun(@(t) t{1}{1}, tok(~cellfun(@isempty,tok)), 'UniformOutput', false);
            app.avail_pts = bases;
            app.sel_order = {};
            app.selList.Items = {};
            buildCheckboxes(bases);
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
        existing = cellfun(@entry_label, app.sel_order, 'UniformOutput', false);
        for k = 1:numel(sel)
            if ~any(strcmp(entry_label(sel{k}), existing))
                app.sel_order{end+1} = sel{k};
                existing{end+1}      = entry_label(sel{k}); %#ok<AGROW>
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
        existing = cellfun(@entry_label, app.sel_order, 'UniformOutput', false);
        if ~any(strcmp(lbl, existing))
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
        n = numel(app.sel_order);
        if n == 0
            app.selList.Items = {};
            return;
        end
        lbls = cellfun(@entry_label, app.sel_order, 'UniformOutput', false);
        % Tag first and last with role indicators
        if n == 1
            lbls{1} = ['[HEAD=TAIL] ' lbls{1}];
        else
            lbls{1}   = ['[HEAD] ' lbls{1}];
            lbls{end} = ['[TAIL] ' lbls{end}];
        end
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
            switch app.fmt
                case 'DLC'
                    app.fp = load_fish_points(csvPath);

                case 'named'
                    if numel(app.sel_order) < 3
                        uialert(fig, 'Select at least 3 points (head, middle(s), tail).', 'Points'); return
                    end
                    app.fp = load_fish_points_named(csvPath, app.sel_order, []);
                    app.fp = app.fp(1);

                case 'numbered'
                    % Format C: auto-loads all numbered points via load_fish_points
                    app.fp = load_fish_points(csvPath);

                case 'curves'
                    % Format E: CURVES — one file per fish, may be multi-selected
                    paths = strsplit(csvPath, ';');
                    fp_arr = load_fish_curves(strtrim(paths{1}));
                    for k_idx = 2:numel(paths)
                        fp_arr(k_idx) = load_fish_curves(strtrim(paths{k_idx})); %#ok<AGROW>
                    end
                    app.fp = fp_arr;

                case 'dual_camera'
                    % Format D: load via load_fish_points_named, which auto-detects dual-cam
                    fp_tmp = load_fish_points_named(csvPath, [], []);
                    app.fp = fp_tmp;
                    app.pect_data = fp_tmp.pect_phase_result;
                    setStatus('Dual-camera loaded. Displaying pectoral phase...');
                    displayPectPhase(app.pect_data);
                    % For kinematics, use cam1 lateral points (already in fp.points)
                    % Fall through to transform + kinematics below

                otherwise
                    uialert(fig, 'Unknown format — browse a file first.', 'Format'); return
            end
        catch ME; uialert(fig, ME.message, 'Load error'); setStatus('Load failed.'); return; end

        setStatus('Remapping axes...');
        try
            app.fp = remap_axes(app.fp);
        catch ME; uialert(fig, ME.message, 'Axis remap error'); setStatus('Failed.'); return; end

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

        % Stage all animals from this run for CSV export
        collect_kine_rows();
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

        % Identify which tracked points are head and tail
        if isfield(fp, 'point_names') && ~isempty(fp.point_names)
            pnames   = fp.point_names;
            headName = pnames{1};
            tailName = pnames{end};
            nPts     = numel(pnames);
            L{end+1} = sprintf('HEAD point (pt 1/%d):  %s', nPts, headName);
            L{end+1} = sprintf('TAIL point (pt %d/%d): %s', nPts, nPts, tailName);
        end

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
        % Body axis reference line at Y=0 (fish is stationary, body anchored)
        yline(app.ax(1), 0, 'k--', 'LineWidth', 1.2, 'Alpha', 0.5, ...
              'DisplayName', 'Body axis (Y = 0)');
        % Head anchor marker at origin
        xline(app.ax(1), 0, 'Color', [0.6 0.6 0.6], 'LineWidth', 0.8, ...
              'Alpha', 0.4, 'HandleVisibility', 'off');

        for ii = 1:numel(idx_show)
            f = valid(idx_show(ii));
            plot(app.ax(1), k.X_interp(f,:), k.Y_interp(f,:), ...
                 'Color', [cmap(ii,:) alpha_val], 'LineWidth', lw);
        end

        % Raw measured points at middle valid frame
        mf = valid(round(end/2));
        scatter(app.ax(1), fp.X(mf,:), fp.Y(mf,:), 50, 'r', 'filled', ...
                'DisplayName', 'Tracked points (mid frame)', 'HandleVisibility', 'on');

        % Colorbar: frame number maps to time
        colormap(app.ax(1), cmap);
        cb = colorbar(app.ax(1));
        cb.Label.String = 'Frame  (blue = early  →  yellow = late)';
        cb.FontSize = 8;
        clim(app.ax(1), [valid(idx_show(1)), valid(idx_show(end))]);

        % X axis: fixed 0 (head) to 1 (tail) in BL — fish is stationary so
        % body length is normalized. Slight padding reveals head/tail markers.
        xlim(app.ax(1), [-0.05, 1.10]);

        % Y axis: symmetric about 0, minimum ±1 BL, expands if data exceeds that.
        % This ensures consistent scale across formats and animals.
        all_y = k.Y_interp(~isnan(k.Y_interp));
        if ~isempty(all_y)
            y_extent = max(abs(all_y));
            y_lim    = max(y_extent * 1.20, 1.0);   % at least ±1 BL, 20% pad beyond data
        else
            y_lim = 1.0;
        end
        ylim(app.ax(1), [-y_lim, y_lim]);

        hold(app.ax(1),'off');
        xlabel(app.ax(1), 'Body position  (BL,  0 = head  →  1 = tail)');
        ylabel(app.ax(1), 'Lateral displacement  (BL,  0 = body axis)');

        % Annotate head and tail point names directly on the plot
        if isfield(fp, 'point_names') && ~isempty(fp.point_names)
            pnames = fp.point_names;
            text(app.ax(1), 0, 0, sprintf('  HEAD\n  (%s)', pnames{1}), ...
                 'FontSize', 8, 'Color', [0.15 0.45 0.15], 'FontWeight', 'bold', ...
                 'VerticalAlignment', 'bottom', 'Clipping', 'on');
            text(app.ax(1), 1, 0, sprintf('TAIL  \n(%s)  ', pnames{end}), ...
                 'FontSize', 8, 'Color', [0.6 0.15 0.15], 'FontWeight', 'bold', ...
                 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', 'Clipping', 'on');
        end

        nSkipped = sum(any(isnan(k.Y_interp), 2));
        title(app.ax(1), ...
            sprintf(['Body-axis midlines  —  fish stationary, swimming against flow\n' ...
                     'n = %d frames shown, %d skipped (occluded)  |  ' ...
                     'color = time progression'], ...
                    numel(idx_show), nSkipped), ...
            'FontSize', 9);
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

        % Allow 3D named (Format B) or 3D numbered (Format C)
        is_3d_ok = isfield(app.fp, 'has_z') && app.fp(1).has_z;
        if ~is_3d_ok
            uialert(fig, '3D data required for fin analysis. Load a file with X, Y, Z columns.', '3D required');
            return
        end

        rootName = app.finRootDrop.Value;
        tipName  = app.finTipDrop.Value;

        if strcmp(rootName, tipName)
            uialert(fig, 'Root and tip must be different points.', 'Selection'); return
        end
        if contains(rootName, 'N/A') || contains(tipName, 'N/A')
            uialert(fig, 'Named or numbered 3D format required for fin analysis.', 'Format'); return
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

        % Stage fin row for CSV export
        collect_fin_row(fd);

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
    %  PECTORAL PHASE DISPLAY  (Format D dual-camera)
    %% ================================================================

    function displayPectPhase(pr)
        if isempty(pr) || ~isstruct(pr)
            app.pectResultArea.Value = {'No pectoral phase data. Load a dual-camera (Format D) CSV.'};
            return;
        end

        % ---- Results text ----
        R = {};
        R{end+1} = '=== PECTORAL FIN PHASE ANALYSIS ===';
        R{end+1} = sprintf('Right fin (pt2) camera:   %s', pr.cam_right);
        R{end+1} = sprintf('Left  fin (pt12) camera:  %s', pr.cam_left);
        R{end+1} = repmat('-', 1, 50);
        if isnan(pr.phase_shift_deg)
            R{end+1} = sprintf('Result: %s', pr.classification);
        else
            R{end+1} = sprintf('CLASSIFICATION:  >>> %s <<<', upper(pr.classification));
            R{end+1} = repmat('-', 1, 50);
            R{end+1} = sprintf('Phase shift (right -> left):  %.1f deg', pr.phase_shift_deg);
            R{end+1} = sprintf('Phase right fin:              %.1f deg', pr.phase_right_deg);
            R{end+1} = sprintf('Phase left  fin:              %.1f deg', pr.phase_left_deg);
            R{end+1} = sprintf('Peak cross-corr lag:          %+d frames', pr.peak_lag_frames);
            R{end+1} = sprintf('Dominant freq (normalized):   %.4f cycles/frame', pr.dom_freq_norm);
            R{end+1} = sprintf('Valid frames:                 %d', pr.n_valid);
            R{end+1} = repmat('-', 1, 50);
            R{end+1} = '  In-phase    : |shift| <= 45 deg or >= 315 deg';
            R{end+1} = '  Antiphase   : 135 deg <= shift <= 225 deg';
            R{end+1} = '  Intermediate: all other values';
        end
        app.pectResultArea.Value = R;

        if isnan(pr.phase_shift_deg) || ~isfield(pr, 'sig_right')
            return;
        end

        frames = 1:numel(pr.sig_right);

        % ---- Plot 1: Time-domain signals ----
        cla(app.pectAx(1)); hold(app.pectAx(1), 'on');
        plot(app.pectAx(1), frames, pr.sig_right, 'b-', 'LineWidth', 1.4, ...
             'DisplayName', sprintf('Right pect tip (pt2) [%s]', pr.cam_right));
        plot(app.pectAx(1), frames, pr.sig_left,  'r-', 'LineWidth', 1.4, ...
             'DisplayName', sprintf('Left  pect base (pt12) [%s]', pr.cam_left));
        yline(app.pectAx(1), 0, 'k:', 'LineWidth', 0.8, 'HandleVisibility', 'off');
        legend(app.pectAx(1), 'Location', 'best', 'FontSize', 9);
        title(app.pectAx(1), sprintf('Pectoral fin signals  |  Classification: %s  |  Phase shift: %.1f deg', ...
              pr.classification, pr.phase_shift_deg), 'FontSize', 10);
        xlabel(app.pectAx(1), 'Frame');
        ylabel(app.pectAx(1), 'Y position (pixels)');
        hold(app.pectAx(1), 'off');

        % ---- Plot 2: Cross-correlation ----
        cla(app.pectAx(2));
        lags = pr.xcorr_lags;
        cc   = pr.xcorr_vals;
        [~, pk_idx] = max(cc);
        hold(app.pectAx(2), 'on');
        fill(app.pectAx(2), [lags fliplr(lags)], [cc zeros(size(cc))], ...
             [0.6 0.8 1], 'EdgeColor', 'none', 'FaceAlpha', 0.5, 'HandleVisibility', 'off');
        plot(app.pectAx(2), lags, cc, 'b-', 'LineWidth', 1.8);
        xline(app.pectAx(2), 0, 'k--', 'LineWidth', 1, 'HandleVisibility', 'off');
        scatter(app.pectAx(2), lags(pk_idx), cc(pk_idx), 80, 'r', 'filled', ...
                'DisplayName', sprintf('Peak lag = %+d frames', pr.peak_lag_frames));
        legend(app.pectAx(2), 'Location', 'best', 'FontSize', 9);
        title(app.pectAx(2), 'Cross-correlation: right vs left pectoral fin', 'FontSize', 10);
        xlabel(app.pectAx(2), 'Lag (frames, positive = left lags right)');
        ylabel(app.pectAx(2), 'Normalized cross-correlation');
        hold(app.pectAx(2), 'off');
    end

    %% ================================================================
    %  AXIS MAPPING HELPERS
    %% ================================================================

    function validateAxisMap()
        tx = app.axTopXDrop.Value;
        ty = app.axTopYDrop.Value;
        tz = app.axVertDrop.Value;
        % Strip signs to get base axes
        bx = tx(2); by = ty(2); bz = tz(2);
        if numel(unique({bx,by,bz})) < 3
            app.axWarnLabel.Text = sprintf('Warning: duplicate axes (%s,%s,%s) — mapping invalid', bx,by,bz);
            app.axWarnLabel.FontColor = [0.75 0.2 0.1];
        else
            app.axWarnLabel.Text = sprintf('OK: CSV %s→topX  %s→topY  %s→vert', tx, ty, tz);
            app.axWarnLabel.FontColor = [0.1 0.55 0.2];
        end
    end

    function fp = remap_axes(fp)
    % REMAP_AXES  Permute and negate fp.points columns so that after remapping:
    %   dim-1 = user-defined top-view X,  dim-2 = top-view Y,  dim-3 = vertical Z.
    % Only acts when the file has 3 dimensions AND the mapping is not the default identity.

        tx = app.axTopXDrop.Value;   % e.g. '+Z'
        ty = app.axTopYDrop.Value;   % e.g. '-X'
        tz = app.axVertDrop.Value;   % e.g. '+Y'

        % Map axis letter to column index in .points dim-3
        ax2col = struct('X',1,'Y',2,'Z',3);

        col_x = ax2col.(tx(2));  sgn_x = 2*strcmp(tx(1),'+')-1;
        col_y = ax2col.(ty(2));  sgn_y = 2*strcmp(ty(1),'+')-1;
        col_z = ax2col.(tz(2));  sgn_z = 2*strcmp(tz(1),'+')-1;

        % Identity check — skip if default (+X,+Y,+Z)
        is_identity = (col_x==1 && sgn_x==1) && (col_y==2 && sgn_y==1) && (col_z==3 && sgn_z==1);

        for fi = 1:numel(fp)
            if isfield(fp(fi),'pre_transformed') && fp(fi).pre_transformed
                continue;   % CURVES data already in body frame — skip
            end
            nDims = size(fp(fi).points, 3);
            if nDims < 3
                continue;   % 2-D data has no Z to remap
            end
            if is_identity
                continue;
            end
            pts = fp(fi).points;  % [nFrames x nPoints x 3]
            new_pts = NaN(size(pts));
            new_pts(:,:,1) = sgn_x .* pts(:,:,col_x);
            new_pts(:,:,2) = sgn_y .* pts(:,:,col_y);
            new_pts(:,:,3) = sgn_z .* pts(:,:,col_z);
            fp(fi).points = new_pts;
        end
    end

    %% ================================================================
    %  CSV EXPORT HELPERS
    %% ================================================================

    function collect_kine_rows()
    % Build one flat row per animal from the current kine results and stage it.
        if isempty(app.kine) || isempty(app.fp), return; end
        timestamp = datestr(now, 'yyyy-mm-dd HH:MM:SS'); %#ok<TNOW1,DATST>
        srcFile   = strtrim(app.fileField.Value);
        fps_val   = app.fpsField.Value;

        for fi = 1:numel(app.kine)
            k  = app.kine(fi);
            fp = app.fp(fi);
            has_z = isfield(fp,'has_z') && fp.has_z;

            r = struct();
            % ---- Provenance ----
            r.timestamp        = timestamp;
            r.source_file      = srcFile;
            r.animal           = k.name;
            r.data_format      = app.fmt;
            if isfield(fp, 'format') && ~isempty(fp.format)
                r.data_format  = fp.format;
            end
            r.fps              = fps_val;
            r.n_frames         = size(fp.points, 1);

            % ---- Point identity ----
            if isfield(fp,'point_names') && ~isempty(fp.point_names)
                r.head_point = fp.point_names{1};
                r.tail_point = fp.point_names{end};
                r.n_points   = numel(fp.point_names);
            else
                r.head_point = '';
                r.tail_point = '';
                r.n_points   = NaN;
            end

            % ---- Beat frequency ----
            r.head_TBF_Hz      = k.head_TBF;
            r.tail_TBF_Hz      = k.tail_TBF;
            r.head_TBF_Z_Hz    = NaN;
            r.tail_TBF_Z_Hz    = NaN;
            if has_z
                r.head_TBF_Z_Hz = k.headZ_TBF;
                r.tail_TBF_Z_Hz = k.tailZ_TBF;
            end

            % ---- Lateral (Y) amplitude ----
            r.head_amp_Y_BL        = k.headAmp;
            r.tail_amp_Y_BL        = k.tailAmp;
            r.head_tail_amp_ratio  = k.headTailAmpRatio;
            r.min_amp_Y_BL         = k.minAmp;
            r.min_amp_Y_loc_BL     = k.minAmpLoc;
            r.max_amp_Y_BL         = k.maxAmp;
            r.max_amp_Y_loc_BL     = k.maxAmpLoc;

            % ---- Dorso-ventral (Z) amplitude ----
            r.head_amp_Z_BL    = NaN;
            r.tail_amp_Z_BL    = NaN;
            r.min_amp_Z_BL     = NaN;
            r.min_amp_Z_loc_BL = NaN;
            r.max_amp_Z_BL     = NaN;
            r.max_amp_Z_loc_BL = NaN;
            if has_z
                r.head_amp_Z_BL    = k.headAmpZ;
                r.tail_amp_Z_BL    = k.tailAmpZ;
                r.min_amp_Z_BL     = k.minAmpZ;
                r.min_amp_Z_loc_BL = k.minAmpZLoc;
                r.max_amp_Z_BL     = k.maxAmpZ;
                r.max_amp_Z_loc_BL = k.maxAmpZLoc;
            end

            % ---- Propulsive wave ----
            r.wavelength_BL    = k.wavelength;

            % ---- Curvature ----
            r.max_curv_XY      = k.maxCurv;
            r.max_curv_XY_loc  = k.maxCurvLoc;
            r.max_curv_3D      = NaN;
            r.max_curv_3D_loc  = NaN;
            if has_z
                r.max_curv_3D     = k.maxCurv3D;
                r.max_curv_3D_loc = k.maxCurv3DLoc;
            end

            % ---- Row type tag ----
            r.row_type = 'kinematics';

            app.csv_rows{end+1} = r;
        end
        updateCSVLabel();
    end

    function collect_fin_row(fd)
    % Build one flat row from fin kinematics results and stage it.
        if nargin < 1 || isempty(fd), return; end
        timestamp   = datestr(now, 'yyyy-mm-dd HH:MM:SS'); %#ok<TNOW1,DATST>
        animalName  = '';
        if ~isempty(app.fishList.Items) && ~strcmp(app.fishList.Items{1}, '(run analysis first)')
            animalName = app.fishList.Value;
        end

        r = struct();
        r.timestamp       = timestamp;
        r.source_file     = strtrim(app.fileField.Value);
        r.animal          = animalName;
        r.data_format     = 'fin_3D';
        r.fps             = app.fpsField.Value;
        r.n_frames        = fd.n_frames;
        r.head_point      = fd.rootName;
        r.tail_point      = fd.tipName;
        r.n_points        = 2;

        % Fill kinematics columns with NaN so row merges cleanly
        r.head_TBF_Hz     = NaN; r.tail_TBF_Hz     = NaN;
        r.head_TBF_Z_Hz   = NaN; r.tail_TBF_Z_Hz   = NaN;
        r.head_amp_Y_BL   = NaN; r.tail_amp_Y_BL   = NaN;
        r.head_tail_amp_ratio = NaN;
        r.min_amp_Y_BL    = NaN; r.min_amp_Y_loc_BL = NaN;
        r.max_amp_Y_BL    = NaN; r.max_amp_Y_loc_BL = NaN;
        r.head_amp_Z_BL   = NaN; r.tail_amp_Z_BL   = NaN;
        r.min_amp_Z_BL    = NaN; r.min_amp_Z_loc_BL = NaN;
        r.max_amp_Z_BL    = NaN; r.max_amp_Z_loc_BL = NaN;
        r.wavelength_BL   = NaN;
        r.max_curv_XY     = NaN; r.max_curv_XY_loc  = NaN;
        r.max_curv_3D     = NaN; r.max_curv_3D_loc  = NaN;

        % Fin-specific fields
        r.fin_root             = fd.rootName;
        r.fin_tip              = fd.tipName;
        r.fin_mean_length      = fd.mean_length;
        r.fin_std_length       = fd.std_length;
        r.fin_total_tip_dist   = fd.total_dist;
        r.fin_mean_speed       = fd.mean_speed;
        r.fin_std_speed        = fd.std_speed;
        r.fin_peak_speed       = fd.peak_speed;
        r.fin_mean_yaw_deg     = fd.mean_yaw;
        r.fin_std_yaw_deg      = fd.std_yaw;
        r.fin_range_yaw_deg    = fd.range_yaw;
        r.fin_mean_pitch_deg   = fd.mean_pitch;
        r.fin_std_pitch_deg    = fd.std_pitch;
        r.fin_range_pitch_deg  = fd.range_pitch;
        r.fin_mean_roll_deg    = fd.mean_roll;
        r.fin_std_roll_deg     = fd.std_roll;
        r.fin_range_roll_deg   = fd.range_roll;
        r.fin_mean_ang_vel     = fd.mean_ang_vel;
        r.fin_std_ang_vel      = fd.std_ang_vel;
        r.fin_peak_ang_vel     = fd.peak_ang_vel;
        r.fin_n_valid_frames   = fd.n_valid;
        r.fin_pct_valid        = fd.pct_valid;

        r.row_type = 'fin_kinematics';

        app.csv_rows{end+1} = r;
        updateCSVLabel();
    end

    function updateCSVLabel()
        n = numel(app.csv_rows);
        if n == 0
            app.csvRowLabel.Text = '0 rows staged';
            app.csvRowLabel.FontColor = [0.35 0.35 0.35];
        else
            if n == 1, plural = ''; else, plural = 's'; end
            app.csvRowLabel.Text = sprintf('%d row%s staged', n, plural);
            app.csvRowLabel.FontColor = [0.10 0.48 0.22];
        end
    end

    function onClearCSV()
        if isempty(app.csv_rows), return; end
        sel = uiconfirm(fig, ...
            sprintf('Clear all %d staged rows?', numel(app.csv_rows)), ...
            'Clear staged rows', ...
            'Options', {'Clear','Cancel'}, 'DefaultOption', 2, 'CancelOption', 2);
        if strcmp(sel,'Clear')
            app.csv_rows = {};
            updateCSVLabel();
            setStatus('Staged CSV rows cleared.');
        end
    end

    function onExportCSV()
        if isempty(app.csv_rows)
            uialert(fig, 'No rows staged yet. Run an analysis first.', 'Nothing to export');
            return
        end

        % Ask: new file or append?
        if ~isempty(app.csv_path) && isfile(app.csv_path)
            choice = uiconfirm(fig, ...
                sprintf('Append %d new row(s) to:\n%s\n\nor write a new file?', ...
                        numel(app.csv_rows), app.csv_path), ...
                'Export CSV', ...
                'Options', {'Append to existing', 'New file', 'Cancel'}, ...
                'DefaultOption', 1, 'CancelOption', 3);
            if strcmp(choice, 'Cancel'), return; end
            do_append = strcmp(choice, 'Append to existing');
        else
            choice = uiconfirm(fig, ...
                sprintf('Export %d row(s) to a new CSV file?', numel(app.csv_rows)), ...
                'Export CSV', ...
                'Options', {'Choose file', 'Cancel'}, ...
                'DefaultOption', 1, 'CancelOption', 2);
            if strcmp(choice, 'Cancel'), return; end
            do_append = false;
        end

        % If new file, ask for path
        if ~do_append
            [fn, fp_] = uiputfile('*.csv', 'Save kinematics CSV as...', ...
                                   fullfile(fileparts(strtrim(app.fileField.Value)), ...
                                            'kinemetrix_export.csv'));
            if isequal(fn, 0), return; end
            app.csv_path = fullfile(fp_, fn);
        end

        % Write
        try
            write_csv_rows(app.csv_rows, app.csv_path, do_append);
            n = numel(app.csv_rows);
            app.csv_rows = {};   % clear staged rows after successful write
            updateCSVLabel();
            setStatus(sprintf('Exported %d row(s) → %s', n, app.csv_path));
        catch ME
            uialert(fig, ME.message, 'Export failed');
        end
    end

    function write_csv_rows(rows, outpath, do_append)
    % Union all field names across all rows, write header + data.
    % Missing fields in a row are written as empty (for strings) or NaN (for numbers).
        if isempty(rows), return; end

        % Collect all column names in insertion order
        all_cols = {};
        for ri = 1:numel(rows)
            flds = fieldnames(rows{ri});
            for fi2 = 1:numel(flds)
                if ~any(strcmp(flds{fi2}, all_cols))
                    all_cols{end+1} = flds{fi2}; %#ok<AGROW>
                end
            end
        end
        nCols = numel(all_cols);
        nRows = numel(rows);

        % Build cell matrix [nRows x nCols]
        C = cell(nRows, nCols);
        for ri = 1:nRows
            for ci = 1:nCols
                col = all_cols{ci};
                if isfield(rows{ri}, col)
                    val = rows{ri}.(col);
                    if ischar(val) || isstring(val)
                        C{ri,ci} = char(val);
                    elseif isnumeric(val) && isscalar(val)
                        if isnan(val)
                            C{ri,ci} = '';
                        else
                            C{ri,ci} = num2str(val, '%.6g');
                        end
                    else
                        C{ri,ci} = '';
                    end
                else
                    C{ri,ci} = '';
                end
            end
        end

        % If appending, read existing header to align columns
        if do_append && isfile(outpath)
            fid_r = fopen(outpath, 'r');
            existing_header = strsplit(strtrim(fgetl(fid_r)), ',');
            fclose(fid_r);
            % Merge: existing cols first, then any new ones
            merged = existing_header;
            for ci = 1:nCols
                if ~any(strcmp(all_cols{ci}, merged))
                    merged{end+1} = all_cols{ci}; %#ok<AGROW>
                end
            end
            % Re-map C into merged column order
            C2 = cell(nRows, numel(merged));
            for ci = 1:numel(merged)
                src = find(strcmp(merged{ci}, all_cols), 1);
                if ~isempty(src)
                    C2(:,ci) = C(:,src);
                end
            end
            all_cols = merged;
            C = C2;
            fid = fopen(outpath, 'a');
        else
            fid = fopen(outpath, 'w');
            % Write header
            fprintf(fid, '%s\n', strjoin(all_cols, ','));
        end

        if fid == -1
            error('Cannot open file for writing: %s', outpath);
        end

        % Write data rows
        for ri = 1:nRows
            row_strs = C(ri,:);
            % Quote any cell containing a comma
            for ci2 = 1:numel(row_strs)
                if contains(row_strs{ci2}, ',')
                    row_strs{ci2} = ['"' row_strs{ci2} '"'];
                end
            end
            fprintf(fid, '%s\n', strjoin(row_strs, ','));
        end
        fclose(fid);
    end

    function setStatus(msg)
        app.statusLabel.Text = msg; drawnow;
    end

    function y_out = section_label(parent, txt, y_in)
        uilabel(parent, 'Text', txt, 'Position', [12 y_in-18 316 18], ...
                'FontSize', 10, 'FontWeight', 'bold', 'FontColor', [0.5 0.5 0.5]);
        y_out = y_in - 24;
    end

end