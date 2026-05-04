%% demo_kinematics.m
% Full pipeline: load -> transform -> kinematics -> plots
%
% Figures produced:
%   1. Interpolated midlines overlaid (coloured by frame)
%   2. Mean amplitude envelope along the body
%   3. Curvature profile along the body
%   4. Head and tail FFT power spectra (beat frequencies)

clear; clc;

%% ---- USER SETTINGS ----
CSV_FILE = 'data.csv';   % path to your tracking CSV
FPS      = 100;          % frames per second
MIN_FREQ = 0.5;          % minimum plausible beat frequency (Hz)
FISH_IDX = 1;            % which fish to plot (1-based)
% -------------------------

%% 1. Load, transform, compute kinematics
fp   = load_fish_points(CSV_FILE);
fp   = transform_fish(fp);
kine = compute_kinematics(fp, FPS, MIN_FREQ);

fi = FISH_IDX;
k  = kine(fi);

fprintf('\n=== %s ===\n',          k.name);
fprintf('  Head TBF        : %.3f Hz\n',  k.head_TBF);
fprintf('  Tail TBF        : %.3f Hz\n',  k.tail_TBF);
fprintf('  Head amplitude  : %.4f BL\n',  k.headAmp);
fprintf('  Tail amplitude  : %.4f BL\n',  k.tailAmp);
fprintf('  Head/tail ratio : %.4f\n',     k.headTailAmpRatio);
fprintf('  Min amplitude   : %.4f BL at body pos %.2f\n', k.minAmp, k.minAmpLoc);
fprintf('  Max amplitude   : %.4f BL at body pos %.2f\n', k.maxAmp, k.maxAmpLoc);
fprintf('  Wavelength      : %.4f BL\n',  k.wavelength);
fprintf('  Max curvature   : %.4f at body pos %.2f\n',    k.maxCurv, k.maxCurvLoc);

%% 2. Figure 1 — Interpolated midlines overlaid
valid_idx = find(~any(isnan(k.Y_interp), 2));
cmap      = parula(numel(valid_idx));

figure('Name', sprintf('%s — Interpolated midlines', k.name), ...
       'Position', [50 550 700 300]);
hold on; grid on; axis equal;
yline(0, 'k--', 'LineWidth', 0.8);

for i = 1:numel(valid_idx)
    f = valid_idx(i);
    plot(k.X_interp(f,:), k.Y_interp(f,:), ...
         'Color', [cmap(i,:), 0.35], 'LineWidth', 0.8);
end

% Overlay raw measured points for the middle valid frame
mid_f = valid_idx(round(end/2));
scatter(fp(fi).X(mid_f,:), fp(fi).Y(mid_f,:), 60, 'r', 'filled', ...
        'DisplayName', 'Measured pts (mid frame)');

colormap(parula);
cb = colorbar; cb.Label.String = 'Frame (early → late)';
clim([valid_idx(1), valid_idx(end)]);
xlabel('X (nose = 0, BL)'); ylabel('Y (BL)');
title(sprintf('%s — FFT-interpolated midlines (%d frames)', ...
      k.name, numel(valid_idx)));
legend('Location','best');

%% 3. Figure 2 — Amplitude envelope
figure('Name', sprintf('%s — Amplitude envelope', k.name), ...
       'Position', [50 180 600 300]);
hold on; grid on;

fill([k.s_norm, fliplr(k.s_norm)], ...
     [k.amp_mean + k.amp_std, fliplr(k.amp_mean - k.amp_std)], ...
     [0.6 0.8 1], 'EdgeColor','none', 'FaceAlpha', 0.5);
plot(k.s_norm, k.amp_mean, 'b-', 'LineWidth', 2);

xline(k.minAmpLoc, 'g--', 'LineWidth', 1, 'Label','min amp');
xline(k.maxAmpLoc, 'r--', 'LineWidth', 1, 'Label','max amp');

xlabel('Normalised body position (0 = nose, 1 = tail)');
ylabel('Mean |Y| half-amplitude (BL)');
title(sprintf('%s — Amplitude envelope (± 1 SD)', k.name));
legend({'± 1 SD','Mean'}, 'Location','northwest');

%% 4. Figure 3 — Curvature profile
figure('Name', sprintf('%s — Curvature', k.name), ...
       'Position', [660 180 600 300]);
hold on; grid on;

fill([k.s_norm, fliplr(k.s_norm)], ...
     [k.curv_mean + k.curv_std, fliplr(k.curv_mean - k.curv_std)], ...
     [1 0.8 0.7], 'EdgeColor','none', 'FaceAlpha', 0.5);
plot(k.s_norm, k.curv_mean, 'r-', 'LineWidth', 2);

xline(k.maxCurvLoc, 'k--', 'LineWidth', 1, ...
      'Label', sprintf('max curv (%.3f)', k.maxCurv));

xlabel('Normalised body position (0 = nose, 1 = tail)');
ylabel('Mean curvature (1/BL)');
title(sprintf('%s — Curvature profile (± 1 SD)', k.name));
legend({'± 1 SD','Mean'}, 'Location','northwest');

%% 5. Figure 4 — FFT power spectra (head & tail beat frequency)
figure('Name', sprintf('%s — Beat frequencies', k.name), ...
       'Position', [660 550 700 450]);

subplot(2,1,1);
plot(k.head_fft_freq, k.head_fft_power, 'b-', 'LineWidth', 1.2);
hold on;
xline(MIN_FREQ,    'k:', 'LineWidth', 1,   'Label', sprintf('min %.1f Hz', MIN_FREQ));
xline(k.head_TBF,  'r--','LineWidth', 1.5, 'Label', sprintf('%.2f Hz', k.head_TBF), ...
      'LabelVerticalAlignment','bottom');
xlabel('Frequency (Hz)'); ylabel('Power');
title(sprintf('%s — Head (P1) beat frequency', k.name));
grid on;

subplot(2,1,2);
plot(k.tail_fft_freq, k.tail_fft_power, 'r-', 'LineWidth', 1.2);
hold on;
xline(MIN_FREQ,    'k:', 'LineWidth', 1,   'Label', sprintf('min %.1f Hz', MIN_FREQ));
xline(k.tail_TBF,  'b--','LineWidth', 1.5, 'Label', sprintf('%.2f Hz', k.tail_TBF), ...
      'LabelVerticalAlignment','bottom');
xlabel('Frequency (Hz)'); ylabel('Power');
title(sprintf('%s — Tail (P5) beat frequency', k.name));
grid on;

sgtitle(sprintf('%s — Beat frequency spectra', k.name));