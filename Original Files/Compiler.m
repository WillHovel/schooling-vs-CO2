%% COMPILER - Import and compile CurveMapper midline data
% Translated from R by E. Goerig & T. Castro-Santos
% Reads individual CurveMapper .xls files and compiles into one tidy table

clear; clc;

%% Load metadata spreadsheet
Multisp = readtable('mega zebegg spreadsheet.xlsx', 'FileType', 'spreadsheet');
Multisp.clip_id = int32(Multisp.clip_id);

% Build filename for each clip (same pattern as R: trimmed_filename + .mp4_CURVES.xls)
Multisp.fn = strcat(Multisp.trimmed_filename, '.mp4_CURVES.xls');

%% Set path to CurveMapper output files
curvemapper_dir = 'C:/Users/valen/Desktop/ZEBEGG/swimming/Videos/zebegg curvemapper';

% Body location index (always 1-200, one row per digitized point)
location = (1:200)';

% Maximum number of midlines CurveMapper can output
max_midlines = 13;

%% Preallocate output tables
% frames: stores the frame number for each midline in each clip
% midlines: final long-format tidy table
frames_all  = table();
midlines_all = table();

%% Main loop - one iteration per clip
clip_ids = Multisp.clip_id;

for idx = 1:length(clip_ids)
    i  = clip_ids(idx);
    fn = Multisp.fn{Multisp.clip_id == i};
    filepath = fullfile(curvemapper_dir, fn);
    
    % Read the CurveMapper Excel file (no header - col_names=FALSE in R)
    try
        raw = readmatrix(filepath, 'FileType', 'spreadsheet', ...
                         'NumHeaderLines', 0);
    catch
        warning('Could not read file for clip %d: %s', i, filepath);
        continue;
    end
    
    nc = size(raw, 2); % number of columns in this file
    
    % Pad to 26 columns if fewer midlines were digitized
    % CurveMapper outputs 2 cols per midline (x,y) plus possibly frame row
    % 13 midlines = 26 data columns; pad with NaN if fewer
    if nc < 26
        pad = nan(size(raw, 1), 26 - nc);
        raw = [raw, pad];
    end
    
    %% Extract frame numbers (row 1 of the file)
    % CurveMapper stores frame numbers in row 1, columns 1,3,5,...,25
    % (every other column, skipping the y columns)
    frame_cols = 1:2:26; % columns 1,3,5,7,...,25 -> 13 values
    frame_row  = raw(1, frame_cols); % 1x13 vector of frame numbers
    
    % Build a small table: one row per midline for this clip
    for m = 1:max_midlines
        if ~isnan(frame_row(m))
            row = table(int32(i), m, frame_row(m), ...
                'VariableNames', {'clip_id','midline','frame'});
            frames_all = [frames_all; row];
        end
    end
    
    %% Extract x,y coordinate data (rows 2-201)
    % Each pair of columns is one midline: col1=x01, col2=y01, col3=x02...
    data_rows = raw(2:201, :); % 200 rows x 26 cols
    
    for m = 1:max_midlines
        x_col = (m-1)*2 + 1; % column index for x of midline m
        y_col = (m-1)*2 + 2; % column index for y of midline m
        
        x_vals = data_rows(:, x_col); % 200x1
        y_vals = data_rows(:, y_col); % 200x1
        
        % Skip if entire midline is NaN (was padded — not digitized)
        if all(isnan(x_vals))
            continue;
        end
        
        % Build tidy rows for this clip x midline
        n_locs = length(location);
        clip_col = repmat(int32(i), n_locs, 1);
        mid_col  = repmat(m, n_locs, 1);
        
        chunk = table(clip_col, mid_col, location, x_vals, y_vals, ...
            'VariableNames', {'clip_id','midline','location','x','y'});
        
        midlines_all = [midlines_all; chunk];
    end
    
    fprintf('Clip %d compiled (%d rows added)\n', i, max_midlines * 200);
end

%% Remove rows where x or y are NaN
% (these came from padded midlines that weren't actually digitized)
valid = ~isnan(midlines_all.x) & ~isnan(midlines_all.y);
midlines_all = midlines_all(valid, :);

%% Merge frame numbers into midlines table
midlines_all = innerjoin(midlines_all, frames_all, 'Keys', {'clip_id','midline'});

%% Reorder columns to match R output: clip_id, midline, location, frame, x, y
midlines_all = midlines_all(:, {'clip_id','midline','location','frame','x','y'});

%% Sort by clip_id, midline, location
midlines_all = sortrows(midlines_all, {'clip_id','midline','location'});

%% Diagnostics
fprintf('\nTotal rows: %d\n', height(midlines_all));
fprintf('Unique clips: %d\n', length(unique(midlines_all.clip_id)));

% Distribution of number of midlines per clip
nb_midlines = varfun(@max, midlines_all, 'InputVariables', 'midline', ...
    'GroupingVariables', 'clip_id');
figure;
histogram(nb_midlines.max_midline);
xlabel('Number of midlines per clip');
ylabel('Count');
title('Distribution of midline count per clip');

%% Save output
writetable(midlines_all, 'MidlinesZebEgg.csv');
fprintf('\nSaved MidlinesZebEgg.csv\n');