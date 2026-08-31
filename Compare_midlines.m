%% compare_midline_sources.m
% Cross-validates two independent tracking pipelines for the SAME fish/video:
%   Source 1: DLTdv CURVES export   (FF1_..._CURVES.xls)
%   Source 2: DeepLabCut raw export (FF1_..._resnet50....csv)
%
% Produces:
%   1. Overlaid midline plot for each source (matches the style of the
%      reference figure: X = body position 0-1 BL, Y = lateral disp. BL)
%   2. Side-by-side tail-beat frequency / amplitude / wavelength comparison
%
% EDIT THESE THREE LINES before running:
CURVES_FILE   = 'FF1_17_NOV_2025_1.50_S_V_000000_S0001_S0001_avi.avi_CURVES.xls';
DLC_FILE      = 'FF1_17_NOV_2025_1.50_S_V_000000_S0001DLC_resnet50_frogfish_ruby_willJul20shuffle1_100000.csv';
CURVES_FPS    = 30.0;   % <-- FILL IN: fps used for the DLTdv/CURVES digitizing
DLC_FPS       = 30.0;  % confirmed from the DLC _meta.pickle, no need to change
MIN_FREQ      = 0.5;   % Hz, same floor used elsewhere in this toolkit

if isnan(CURVES_FPS)
    error(['Set CURVES_FPS at the top of this script before running: ' ...
           'frequency (Hz) cannot be computed without it.']);
end

%% ---- Load Source 1: CURVES (pre-transformed, skip transform_fish) ----
fp_curves = load_fish_curves(CURVES_FILE);
kine_curves = compute_kinematics(fp_curves, CURVES_FPS, MIN_FREQ);

%% ---- Load Source 2: raw DeepLabCut export ----
% Midline points only (Mid1..Mid7), pectoral points aren't part of the
% body midline, so they're excluded from this comparison. Verify Mid1 is
% truly the head before trusting this order (see prior evaluation notes).
fp_dlc = load_fish_points_named(DLC_FILE, ...
            {'Mid1','Mid2','Mid3','Mid4','Mid5','Mid6','Mid7'}, ...
            [1 2 3 4 5 6 7]);
fp_dlc = transform_fish(fp_dlc);
kine_dlc = compute_kinematics(fp_dlc, DLC_FPS, MIN_FREQ);

%% ---- Plot 1: Overlaid midlines, same style as the reference figure ----
figure('Name', 'Midline comparison', 'Position', [100 100 1200 800]);

subplot(2,1,1);
hold on;
n_show = min(30, size(kine_curves.Y_interp, 1));  % cap for readability
idx = round(linspace(1, size(kine_curves.Y_interp,1), n_show));
for i = idx
    plot(kine_curves.s_norm, kine_curves.Y_interp(i,:), 'LineWidth', 1);
end
title(sprintf('Source 1: DLTdv CURVES  (%.0f fps, %d frames, %d shown)', ...
      CURVES_FPS, size(kine_curves.Y_interp,1), n_show));
xlabel('Body position (BL, 0=head \rightarrow 1=tail)');
ylabel('Lateral displacement (BL)');
ylim([-0.3 0.3]); grid on; box on;

subplot(2,1,2);
hold on;
n_show2 = min(30, size(kine_dlc.Y_interp, 1));
idx2 = round(linspace(1, size(kine_dlc.Y_interp,1), n_show2));
for i = idx2
    plot(kine_dlc.s_norm, kine_dlc.Y_interp(i,:), 'LineWidth', 1);
end
title(sprintf('Source 2: DeepLabCut  (%.0f fps, %d frames, %d shown)', ...
      DLC_FPS, size(kine_dlc.Y_interp,1), n_show2));
xlabel('Body position (BL, 0=head \rightarrow 1=tail)');
ylabel('Lateral displacement (BL)');
ylim([-0.3 0.3]); grid on; box on;

%% ---- Table: frequency / amplitude / wavelength comparison ----
fprintf('\n=== MIDLINE SOURCE COMPARISON: %s ===\n\n', 'FF1_17_NOV_2025');
fprintf('%-28s %14s %14s %10s\n', 'Metric', 'CURVES', 'DeepLabCut', 'Diff %');
fprintf('%s\n', repmat('-', 1, 68));

metrics = {
    'head_TBF (Hz)',        kine_curves.head_TBF,        kine_dlc.head_TBF;
    'tail_TBF (Hz)',        kine_curves.tail_TBF,         kine_dlc.tail_TBF;
    'headAmp (BL)',         kine_curves.headAmp,          kine_dlc.headAmp;
    'tailAmp (BL)',         kine_curves.tailAmp,          kine_dlc.tailAmp;
    'headTailAmpRatio',     kine_curves.headTailAmpRatio, kine_dlc.headTailAmpRatio;
    'wavelength (BL)',      kine_curves.wavelength,       kine_dlc.wavelength;
    'maxCurv (1/BL)',       kine_curves.maxCurv,          kine_dlc.maxCurv;
};

for i = 1:size(metrics,1)
    name = metrics{i,1};  v1 = metrics{i,2};  v2 = metrics{i,3};
    diff_pct = 100 * (v2 - v1) / v1;
    fprintf('%-28s %14.4f %14.4f %9.1f%%\n', name, v1, v2, diff_pct);
end

fprintf('\nNote: frame counts differ (CURVES=%d frames, DLC=%d frames), this\n', ...
        size(kine_curves.Y_interp,1), size(kine_dlc.Y_interp,1));
fprintf('affects duration/statistical power but not the Hz frequency values\n');
fprintf('themselves, as long as CURVES_FPS above is correct.\n');