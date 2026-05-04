%% demo_transform.m
% Loads fish tracking data, applies transform_fish(), and plots a before/after
% comparison for one fish across several frames.

%% 1. Load & transform
fp = load_fish_points('C:\Users\willh\Desktop\Fish_analysis_v2\data.csv');   % <-- change to your CSV path
fp = transform_fish(fp);

%% 2. Pick a fish and a handful of frames to overlay
fish_idx  = 1;                        % which fish to plot (1-based)
n_frames  = 10;                       % how many frames to overlay
frame_step = max(1, floor(size(fp(fish_idx).points, 1) / n_frames));
frame_ids  = 1:frame_step:size(fp(fish_idx).points, 1);
frame_ids  = frame_ids(1:min(n_frames, end));

cmap = parula(numel(frame_ids));      % one colour per frame

%% 3. Plot: before (raw) on the left, after (rotated+translated) on the right
figure('Name', sprintf('Transform check — %s', fp(fish_idx).name), ...
       'Position', [100 100 1000 420]);

% ---------- Before ----------
subplot(1, 2, 1);
hold on; grid on; axis equal;
title(sprintf('%s — Raw points', fp(fish_idx).name));
xlabel('x (pixels)');  ylabel('y (pixels)');

for k = 1:numel(frame_ids)
    f  = frame_ids(k);
    xr = fp(fish_idx).points(f, :, 1);
    yr = fp(fish_idx).points(f, :, 2);
    if any(isnan(xr) | isnan(yr)), continue; end

    plot(xr, yr, '-o', 'Color', cmap(k,:), ...
         'MarkerFaceColor', cmap(k,:), 'MarkerSize', 5, 'LineWidth', 1.2);
    % Mark head (P1) and tail (P5)
    plot(xr(1),   yr(1),   'k^', 'MarkerFaceColor','k',   'MarkerSize', 7); % head
    plot(xr(end), yr(end), 'ks', 'MarkerFaceColor','none','MarkerSize', 7); % tail
end
legend({'midline','','head (P1)','tail (P5)'}, 'Location','best');

% ---------- After ----------
subplot(1, 2, 2);
hold on; grid on; axis equal;
title(sprintf('%s — Rotated & translated (X, Y)', fp(fish_idx).name));
xlabel('X (nose = 0)');  ylabel('Y');
yline(0, 'k--', 'LineWidth', 0.8);   % x-axis reference

for k = 1:numel(frame_ids)
    f  = frame_ids(k);
    xr = fp(fish_idx).X(f, :);
    yr = fp(fish_idx).Y(f, :);
    if any(isnan(xr) | isnan(yr)), continue; end

    plot(xr, yr, '-o', 'Color', cmap(k,:), ...
         'MarkerFaceColor', cmap(k,:), 'MarkerSize', 5, 'LineWidth', 1.2);
    plot(xr(1),   yr(1),   'k^', 'MarkerFaceColor','k',   'MarkerSize', 7);
    plot(xr(end), yr(end), 'ks', 'MarkerFaceColor','none','MarkerSize', 7);
end

colormap(parula);
cb = colorbar;
cb.Label.String = 'Frame (early → late)';
clim([frame_ids(1) frame_ids(end)]);

sgtitle('Fish midline transformation — middle 3 points used for rotation');