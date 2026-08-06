%% Rotating and translating midlines for kinematics analysis
% Translated from R code by T. Castro-Santos & E. Goerig (Sept 23, 2017)
% Adjusts for camera rotation and centers nose position at (0,0)

clear; clc;

%% Load data
midlines = readtable('MidlinesZebEgg.csv');
midlines = midlines(:, 2:7); % Select columns 2-7 (same as R's select)

Multisp = readtable('mega zebegg spreadsheet.xlsx', 'FileType', 'spreadsheet');
Multisp.clip_id = int32(Multisp.clip_id);

%% STANDARDIZATION
% Extract relevant columns from Multisp
TL = Multisp(:, {'clip_id', 'TL', 'species'});
TL.TL = str2double(string(TL.TL)); % Ensure numeric

% Replace NaN TL with 1
TL.TL(isnan(TL.TL)) = 1;

% Create TL2: force TL=1 for specific clip IDs (already digitized at 1:1 scale)
specialClips = [8:11, 16:18, 30, 49, 62, 101:121, 144:146, 155:157, 160:187];
TL.TL2 = TL.TL; % copy TL into TL2
TL.TL2(ismember(TL.clip_id, specialClips)) = 1;

% Join TL info into midlines table (equivalent to inner_join by clip_id)
midlines = innerjoin(midlines, TL, 'Keys', 'clip_id');
midlines.TL = midlines.TL2; % Use TL2 as the working TL

% Fix data entry error in clip 11 (coordinates were 100x too large)
clip11 = midlines.clip_id == 11;
midlines.x(clip11) = midlines.x(clip11) / 100;
midlines.y(clip11) = midlines.y(clip11) / 100;

% Standardize x and y to body lengths (BL)
midlines.xBL = midlines.x ./ midlines.TL2;
midlines.yBL = midlines.y ./ midlines.TL2;

%% ROTATION

% Fix clip 116 (Tuna) - was digitized in reverse, so mirror x-coordinates
clip116 = midlines.clip_id == 116;
midlines.xBL(clip116) = -midlines.xBL(clip116) + 0.5;

midlines.clip_id = double(midlines.clip_id);

% Select subset of midline points: location between 10 and 190
% (roughly the lead 90% of body, excluding tail tip)
chop = midlines(midlines.location >= 10 & midlines.location <= 190, :);
clipid = unique(chop.clip_id);

% Preallocate models table
models = table('Size', [0, 4], ...
    'VariableTypes', {'double','double','double','double'}, ...
    'VariableNames', {'clip_id','a','b','R2'});

%% Fit linear regression for each clip to find swimming angle
for i = 1:length(clipid)
    id = clipid(i);
    a_data = chop(chop.clip_id == id, :);
    
    % Fit linear model: yBL ~ xBL
    coeffs = polyfit(a_data.xBL, a_data.yBL, 1); % returns [slope, intercept]
    b = coeffs(1);       % slope
    a = coeffs(2);       % intercept
    
    % Calculate R-squared manually
    y_pred = polyval(coeffs, a_data.xBL);
    ss_res = sum((a_data.yBL - y_pred).^2);
    ss_tot = sum((a_data.yBL - mean(a_data.yBL)).^2);
    R2 = 1 - ss_res / ss_tot;
    
    models = [models; {id, a, b, R2}];
end

%% Calculate rotation angles
models.alpha = atan(models.b);           % angle from slope
models.theta = 2*pi - models.alpha;      % adjust direction (clockwise rotation)

%% Apply rotation to all midline points
% Rotation formulas (rotating about the y-intercept, not the origin):
%   x' = x*cos(theta) - (y - a)*sin(theta)
%   y' = x*sin(theta) + (y - a)*cos(theta) + a

midlines.xrot = nan(height(midlines), 1);
midlines.yrot = nan(height(midlines), 1);

for i = 1:length(clipid)
    id = clipid(i);
    a     = models.a(models.clip_id == id);
    theta = models.theta(models.clip_id == id);
    
    mask = midlines.clip_id == id;
    x = midlines.xBL(mask);
    y = midlines.yBL(mask);
    
    midlines.xrot(mask) = x .* cos(theta) - (y - a) .* sin(theta);
    midlines.yrot(mask) = x .* sin(theta) + (y - a) .* cos(theta) + a;
end

%% Save rotated plots as PNG files
mkdir('PlotsRot'); % create output folder if it doesn't exist

for i = 1:length(clipid)
    id = clipid(i);
    a_data = midlines(midlines.clip_id == id, :);
    
    fig = figure('Visible', 'off'); % don't display, just save
    hold on;
    
    % Plot each midline as a separate colored line (equivalent to geom_path)
    midline_ids = unique(a_data.midline);
    colors = lines(length(midline_ids)); % colormap with enough colors
    for m = 1:length(midline_ids)
        mid_data = a_data(a_data.midline == midline_ids(m), :);
        plot(mid_data.xrot, mid_data.yrot, 'Color', colors(m, :));
    end
    
    xlim([-0.5, 1.5]);
    ylim([-1, 1]);
    title(['Clip ID - ', num2str(id)]);
    set(gca, 'Color', 'none'); % transparent background equivalent
    box off;
    
    saveas(fig, fullfile('PlotsRot', [num2str(id), '.png']));
    close(fig);
end


%% TRANSLATION
% Shift midlines so the nose (minimum x after rotation) is at x = 0

mkdir('PlotsAdj'); % create output folder if needed

% For each clip+midline combo, find the minimum xrot (nose position)
% then subtract it to translate so nose is at x = 0

midlines.X = nan(height(midlines), 1);
midlines.Y = midlines.yrot; % Y is unchanged — fish oscillate freely in y

clip_midline_combos = unique(midlines(:, {'clip_id', 'midline'}), 'rows');

for i = 1:height(clip_midline_combos)
    id  = clip_midline_combos.clip_id(i);
    mid = clip_midline_combos.midline(i);
    
    mask = midlines.clip_id == id & midlines.midline == mid;
    minx = min(midlines.xrot(mask));
    midlines.X(mask) = midlines.xrot(mask) - minx;
end

%% Save translated plots as PNG
for i = 1:length(clipid)
    id = clipid(i);
    a_data = midlines(midlines.clip_id == id, :);
    
    fig = figure('Visible', 'off');
    hold on;
    
    midline_ids = unique(a_data.midline);
    colors = lines(length(midline_ids));
    for m = 1:length(midline_ids)
        mid_data = a_data(a_data.midline == midline_ids(m), :);
        plot(mid_data.X, mid_data.Y, 'Color', colors(m, :));
    end
    
    xlim([-0.5, 1.5]);
    ylim([-1, 1]);
    title(['Clip ID - ', num2str(id)]);
    set(gca, 'Color', 'none');
    box off;
    
    saveas(fig, fullfile('PlotsAdj', [num2str(id), '.png']));
    close(fig);
end

%% Save final output to CSV
writetable(midlines, 'MidlinesZebEggAdj.csv');

% END



%% PRE-ROTATION DIAGNOSTIC PLOTS -- TESTING ONLY
% Plots raw (xBL, yBL) before rotation — compare against PlotsRot to verify correction

mkdir('Plots'); % create output folder if needed

for i = 1:length(clipid)
    id = clipid(i);
    a_data = chop(chop.clip_id == id, :); % use 'chop' (location 10-190), same as R

    fig = figure('Visible', 'off');
    hold on;

    midline_ids = unique(a_data.midline);
    colors = lines(length(midline_ids));
    for m = 1:length(midline_ids)
        mid_data = a_data(a_data.midline == midline_ids(m), :);
        plot(mid_data.xBL, mid_data.yBL, 'Color', colors(m, :)); % raw BL coords
    end

    xlim([-0.5, 1.5]);
    ylim([-1, 1]);
    title(['Clip ID - ', num2str(id)]);
    set(gca, 'Color', 'none');
    box off;

    saveas(fig, fullfile('Plots', [num2str(id), '.png']));
    close(fig);
end