%% KINEMATICS - Full Analysis
% Translated from R by E. Goerig & T. Castro-Santos (19/07/2021)

clear; clc;

%% Setup & Import data
midlines = readtable('MidlinesAdj.csv');
midlines = midlines(:, 2:17); % Select columns 2-17

Multisp = readtable('Multi-species spreadsheet.xlsx', 'FileType', 'spreadsheet');
Multisp.clip_id = int32(Multisp.clip_id);

% Remove lamprey ammocoetes & transformers
midlines = midlines(~strcmp(midlines.species, 'Petromyzon_marinus_ammocoete'), :);
midlines = midlines(~strcmp(midlines.species, 'Petromyzon_marinus_transformer'), :);
Multisp  = Multisp(~strcmp(Multisp.species,  'Petromyzon_marinus_ammocoete'), :);
Multisp  = Multisp(~strcmp(Multisp.species,  'Petromyzon_marinus_transformer'), :);

% Fix clip 2: midlines 10 and 11 were reversed during digitizing
mask10 = midlines.clip_id == 2 & midlines.midline == 10;
mask11 = midlines.clip_id == 2 & midlines.midline == 11;
midlines.midline(mask10) = 11;
midlines.midline(mask11) = 10;

% Harmonize species names (remove spaces, replace with underscores)
midlines.species = strrep(midlines.species, 'Anguilla rostrata',    'Anguilla_rostrata');
midlines.species = strrep(midlines.species, 'Thunnus albacares',    'Thunnus_albacares');
midlines.species = strrep(midlines.species, 'Caranx hippos',        'Caranx_hippos');
midlines.species = strrep(midlines.species, 'Devario aequipinnatus','Devario_aequipinnatus');

% Tailbeat frequency = fps / totalframes
Multisp.TBF = Multisp.fps ./ Multisp.totalframes;

% Build megakine summary table
megakine = Multisp(:, {'species','clip_id','lab','speedflowBL','speedflow_m','fps','TL','TBF'});
megakine.speedflowBL = double(megakine.speedflowBL);
megakine.TL          = double(megakine.TL);

% Load and merge swimming modes
Swimmingmode = readtable('Swimmingmode.xlsx', 'FileType', 'spreadsheet');
Swimmingmode.Swimmingmode = categorical(Swimmingmode.Swimmingmode);
megakine = innerjoin(megakine, Swimmingmode, 'Keys', 'clip_id');


%% Ground speed & Swim speed

% Get tail tip (max location) at first and last frame for each clip
clipid = unique(midlines.clip_id);

UgBL_vals = nan(length(clipid), 1);

for i = 1:length(clipid)
    id = clipid(i);
    clip_data = midlines(midlines.clip_id == id, :);
    
    min_frame = min(clip_data.frame);
    max_frame = max(clip_data.frame);
    max_loc   = max(clip_data.location);
    
    % Tail tip row at first frame
    a1 = clip_data(clip_data.frame == min_frame & clip_data.location == max_loc, :);
    % Tail tip row at last frame
    a2 = clip_data(clip_data.frame == max_frame & clip_data.location == max_loc, :);
    
    if isempty(a1) || isempty(a2), continue; end
    
    % Distance in BL (absolute x displacement of tail tip)
    dist = abs(a2.xrot(1) - a1.xrot(1));
    
    % Time elapsed in seconds
    fps_i = Multisp.fps(Multisp.clip_id == id);
    time  = (a2.frame(1) - a1.frame(1)) / fps_i;
    
    UgBL_vals(i) = dist / time; % ground speed in BL/s
end

megakine.UgBL = UgBL_vals;
Multisp.UgBL  = UgBL_vals;

% Swim speed = ground speed + flow speed (BL/s)
megakine.UsBL = megakine.UgBL + megakine.speedflowBL;
Multisp.UsBL  = megakine.UsBL;

% Quick diagnostic plots
figure; histogram(megakine.UsBL);
title('Swim speed distribution'); xlabel('UsBL (BL/s)');

figure; 
labs = unique(megakine.lab);
boxplot(megakine.UsBL, megakine.lab);
title('Swim speed by lab'); xlabel('Lab'); ylabel('UsBL (BL/s)');

% Compare video-based swim speed with PIT-based speed (Conte only)
diff_speed = megakine.UsBL - Multisp.UsBL_PIT;
figure; histogram(diff_speed);
title('Difference: video UsBL vs PIT UsBL');

% Save updated Multisp to Excel
writetable(Multisp, 'Multi-Species spreadsheet1.xlsx');


%% Stride Length
% SL = swim speed (BL/s) / tailbeat frequency (TBF)
megakine.SL = megakine.UsBL ./ megakine.TBF;

%% Amplitude
% For each species x clip x body location, find min and max Y across all midlines
% Amplitude = max(Y) - min(Y), HalfAmp = amplitude/2

mkdir('PlotsAmplitude');

% Get unique combinations of species, clip_id, location
[grp, sp_u, id_u, loc_u] = findgroups(...
    midlines.species, midlines.clip_id, midlines.location);

miny = splitapply(@min, midlines.Y, grp);
maxy = splitapply(@max, midlines.Y, grp);

amplitude = table(id_u, sp_u, loc_u, miny, maxy, ...
    'VariableNames', {'clip_id','species','location','miny','maxy'});
amplitude.amplitude = amplitude.maxy - amplitude.miny;
amplitude.HalfAmp   = amplitude.amplitude / 2;
amplitude.location  = double(amplitude.location);

% Save amplitude plots per clip
clipid = unique(midlines.clip_id);

for i = 1:length(clipid)
    id     = clipid(i);
    a_data = amplitude(amplitude.clip_id == id, :);
    if isempty(a_data), continue; end
    
    sp_name = a_data.species{1};
    fig = figure('Visible','off');
    plot(a_data.location, a_data.amplitude, 'k-');
    xlim([-5, 205]); ylim([0, 0.5]);
    xlabel('Location'); ylabel('Amplitude (BL)');
    title([sp_name, ' Clip ID - ', num2str(id)]);
    set(gca, 'Color','none'); box off;
    saveas(fig, fullfile('PlotsAmplitude', [num2str(id), '.png']));
    saveas(fig, fullfile('PlotsAmplitude', [num2str(id), '.pdf']));
    close(fig);
end

%% Head and tail amplitude
head_amp_data = amplitude(amplitude.location >= 1 & amplitude.location <= 5, :);
tail_amp_data = amplitude(amplitude.location > 195, :);

head_ids = unique(head_amp_data.clip_id);
headAmp  = table(head_ids, nan(length(head_ids),1), 'VariableNames',{'clip_id','headAmp'});
for i = 1:length(head_ids)
    id = head_ids(i);
    headAmp.headAmp(i) = mean(head_amp_data.amplitude(head_amp_data.clip_id == id));
end

tail_ids = unique(tail_amp_data.clip_id);
tailAmp  = table(tail_ids, nan(length(tail_ids),1), 'VariableNames',{'clip_id','tailAmp'});
for i = 1:length(tail_ids)
    id = tail_ids(i);
    tailAmp.tailAmp(i) = mean(tail_amp_data.amplitude(tail_amp_data.clip_id == id));
end

megakine.headAmp    = headAmp.headAmp;
megakine.tailAmp    = tailAmp.tailAmp;
megakine.headtailamp = megakine.headAmp ./ megakine.tailAmp;

%% Min and max amplitude location along the body
[~, min_idx] = splitapply(@(x,loc) deal(min(x), loc(x==min(x),1)), ...
    amplitude.amplitude, amplitude.location, findgroups(amplitude.clip_id));
% Simpler explicit loop:
amp_clipids = unique(amplitude.clip_id);
minAmploc = nan(length(amp_clipids),1);
minAmp    = nan(length(amp_clipids),1);
maxAmploc = nan(length(amp_clipids),1);
maxAmp    = nan(length(amp_clipids),1);

for i = 1:length(amp_clipids)
    id     = amp_clipids(i);
    a_data = amplitude(amplitude.clip_id == id, :);
    [minAmp(i),    mi] = min(a_data.amplitude);
    [maxAmp(i),    ma] = max(a_data.amplitude);
    minAmploc(i) = a_data.location(mi);
    maxAmploc(i) = a_data.location(ma);
end

megakine.minAmploc = minAmploc;
megakine.minAmp    = minAmp;
megakine.maxAmploc = maxAmploc;
megakine.maxAmp    = maxAmp;

% Normalized body location (0 to 1)
amplitude.location2 = amplitude.location / 200;

% Save amplitude to CSV
writetable(amplitude, 'Amplitude.csv');

%% Propulsive Wave
% Track peaks moving down the body to calculate wavespeed and wavelength

mkdir('midlineplots');
mkdir('wavespeed');
mkdir('WavespeedAdj');

%% Step 1 - Plot midlines centered on Y for visual inspection
clipid = unique(midlines.clip_id);

for i = 1:length(clipid)
    id     = clipid(i);
    a_data = midlines(midlines.clip_id == id, :);
    
    % Center Y for visualization only
    a_data.Ycenter = a_data.Y - mean(a_data.Y);
    sp_name = a_data.species{1};
    
    mid_ids = unique(a_data.midline);
    n_mids  = length(mid_ids);
    
    fig = figure('Visible','off','Position',[0 0 1400 500]);
    for m = 1:n_mids
        subplot(1, n_mids, m);
        md = a_data(a_data.midline == mid_ids(m), :);
        plot(md.X, md.Ycenter, 'k-');
        xlim([0,1]); ylim([-0.3, 0.3]);
        title(num2str(mid_ids(m)));
    end
    sgtitle([sp_name, ' Clip ID - ', num2str(id)]);
    saveas(fig, fullfile('midlineplots', [num2str(id), '.png']));
    close(fig);
end

%% Step 2 - Identify local peaks (Y maxima) along each midline
% A point is a local peak if Y > both its neighbor before and after

df = midlines;
df.localpeak = zeros(height(df), 1);

clip_mid_combos = unique(df(:, {'clip_id','midline'}), 'rows');

for i = 1:height(clip_mid_combos)
    id  = clip_mid_combos.clip_id(i);
    mid = clip_mid_combos.midline(i);
    mask = df.clip_id == id & df.midline == mid;
    Y_vals = df.Y(mask);
    n = length(Y_vals);
    peaks = zeros(n,1);
    for k = 2:n-1
        dY1 = Y_vals(k) - Y_vals(k-1);
        dY2 = Y_vals(k) - Y_vals(k+1);
        if dY1 > 0 && dY2 > 0
            peaks(k) =  1; % peak
        elseif dY1 < 0 && dY2 < 0
            peaks(k) = -1; % trough
        end
    end
    idx = find(mask);
    df.localpeak(idx) = peaks;
end

YPeak   = df(df.localpeak ==  1, :);
YTrough = df(df.localpeak == -1, :);
YPeak   = YPeak(~any(ismissing(YPeak), 2), :);
YTrough = YTrough(~any(ismissing(YTrough), 2), :);

%% Helper function: track peaks moving anterior to posterior
% (defined below at end of script as a local function)

%% Step 2 - Track peaks down the body (unfiltered)
peakloc = track_peaks(YPeak);

%% Compute Midcount (frame range and proportion of tailbeat)
Midcount = compute_midcount(midlines, Multisp);
peakloc  = innerjoin(peakloc, Midcount, 'Keys', 'clip_id');

% Fix frame offset (clip 115 doesn't start on frame 1)
peakloc.frame   = peakloc.frame - peakloc.minframe + 1;
peakloc.PropTB  = peakloc.frame  ./ peakloc.totalframes;
peakloc.PropLoc = peakloc.location ./ 200;

%% Fit wave regression models (unfiltered) and save plots
wavemodels = fit_wave_models(peakloc, 'wavespeed');

wavemodels = innerjoin(wavemodels, Midcount, 'Keys', 'clip_id');
wavemodels = wavemodels(:, {'clip_id','Intercept','SE_Intercept', ...
    'Slope','SE_Slope','R2','species','nframes','totalframes','Prop_TB'});

% Save unfiltered wave models
writetable(wavemodels, 'TR_Clip_YCenterClip_Nofilter.xlsx');

%% Step 3 - Remove problematic midlines identified by visual QC
a = YPeak;
a = a(~(a.clip_id==1   & a.midline > 7),  :);
a = a(~(a.clip_id==7   & a.midline < 3),  :);
a = a(~(a.clip_id==27  & a.midline == 12), :);
a = a(~(a.clip_id==50  & a.midline == 1),  :);
a = a(~(a.clip_id==65  & a.midline < 4),  :);
a = a(~(a.clip_id==67  & a.midline < 7),  :);
a = a(~(a.clip_id==71  & a.midline < 5),  :);
a = a(~(a.clip_id==72  & a.midline < 4),  :);
a = a(~(a.clip_id==86  & a.midline < 4),  :);
a = a(~(a.clip_id==88  & a.midline < 5),  :);
a = a(~(a.clip_id==130 & a.midline > 7),  :);
a = a(~(a.clip_id==133 & a.midline < 7),  :);
a = a(~(a.clip_id==142 & a.midline < 5),  :);
a = a(~(a.clip_id==154 & a.midline < 11), :);
a = a(~(a.clip_id==234 & a.midline < 3),  :);
a = a(~(a.clip_id==245 & a.midline < 4),  :);
a = a(~(a.clip_id==247 & a.midline == 1),  :);
a = a(~(a.clip_id==249 & a.midline < 4),  :);
a = a(~(a.clip_id==254 & a.midline < 7),  :);

%% Redo peak tracking and wave models on filtered data
peakloc = track_peaks(a);
Midcount = compute_midcount(midlines, Multisp); % recompute from original midlines
peakloc  = innerjoin(peakloc, Midcount, 'Keys', 'clip_id');
peakloc.frame   = peakloc.frame - peakloc.minframe + 1;
peakloc.PropTB  = peakloc.frame  ./ peakloc.totalframes;
peakloc.PropLoc = peakloc.location ./ 200;

wavemodelsAdj = fit_wave_models(peakloc, 'WavespeedAdj');
wavemodelsAdj = innerjoin(wavemodelsAdj, Midcount, 'Keys', 'clip_id');
wavemodelsAdj = wavemodelsAdj(:, {'clip_id','Intercept','SE_Intercept', ...
    'Slope','SE_Slope','R2','species','nframes','totalframes','Prop_TB'});

wavemodelsAdj.wavespeed  = wavemodelsAdj.Slope;       % BL/TB
wavemodelsAdj.wavelength = wavemodelsAdj.wavespeed;   % wavelength = wavespeed * 1 tailbeat

% Merge into megakine
wavelength_tbl = wavemodelsAdj(:, {'clip_id','wavespeed','wavelength'});
megakine = outerjoin(megakine, wavelength_tbl, 'Keys','clip_id','MergeKeys',true);

% Wavespeed in BL/s
megakine.wavespeedBLS = megakine.wavelength .* megakine.TBF;

% Save filtered wave models
writetable(wavemodelsAdj, 'TR_Clip_YCenterClip.xlsx');


%% Reynolds & Strouhal numbers

% Swim speed in m/s (TL in cm -> /100 for meters)
megakine.U = (megakine.UsBL .* megakine.TL) / 100;

% Reynolds number: Re = U * L / v  (v = 1e-6 m2/s for water)
megakine.Re = (megakine.U .* (megakine.TL / 100)) / 1e-6;

% Strouhal number: St = f * A / U
megakine.tailAmp_m = (megakine.tailAmp .* megakine.TL) / 100; % tail amplitude in meters
megakine.St = (megakine.TBF .* megakine.tailAmp_m) ./ megakine.U;

disp(summary(megakine(:, {'Re','St'})));

%% Curvature
% 3-point geometric curvature using lag/lead of 5 locations along body
% Formula uses triangle side lengths A, B, C -> radius of circumscribed circle -> K = 1/R

lagno = 5;
df_curv = midlines;
df_curv.Curv = nan(height(df_curv), 1);

clip_mid_combos = unique(df_curv(:, {'clip_id','midline'}), 'rows');

for i = 1:height(clip_mid_combos)
    id  = clip_mid_combos.clip_id(i);
    mid = clip_mid_combos.midline(i);
    mask = df_curv.clip_id == id & df_curv.midline == mid;
    idx  = find(mask);
    n    = length(idx);
    
    X_vals = df_curv.X(idx);
    Y_vals = df_curv.Y(idx);
    
    for k = lagno+1 : n-lagno
        x1 = X_vals(k - lagno); y1 = Y_vals(k - lagno); % lag point
        x2 = X_vals(k);         y2 = Y_vals(k);          % current
        x3 = X_vals(k + lagno); y3 = Y_vals(k + lagno);  % lead point
        
        A = sqrt((x2-x1)^2 + (y2-y1)^2);
        B = sqrt((x3-x2)^2 + (y3-y2)^2);
        C = sqrt((x3-x1)^2 + (y3-y1)^2);
        D = (A + B + C) / 2;
        
        denom = 4 * sqrt(D*(D-A)*(D-B)*(D-C));
        if denom > 0
            Radius = (A * B * C) / denom;
            df_curv.Curv(idx(k)) = 1 / Radius;
        end
    end
end

% Max curvature per clip (across all midlines and locations)
curv_clipids  = unique(df_curv.clip_id);
maxcurv_vals  = nan(length(curv_clipids), 1);
maxcurv_locs  = nan(length(curv_clipids), 1);
maxcurv_mids  = nan(length(curv_clipids), 1);
maxradius_vals = nan(length(curv_clipids), 1);

for i = 1:length(curv_clipids)
    id     = curv_clipids(i);
    a_data = df_curv(df_curv.clip_id == id, :);
    [maxcurv_vals(i), mi] = max(a_data.Curv);
    maxcurv_locs(i)  = a_data.location(mi);
    maxcurv_mids(i)  = a_data.midline(mi);
    maxradius_vals(i) = 1 / maxcurv_vals(i);
end

maxcurv_tbl = table(curv_clipids, maxcurv_mids, maxcurv_locs, ...
    maxradius_vals, maxcurv_vals, ...
    'VariableNames',{'clip_id','maxcurvmidline','maxcurvloc','maxradius','maxcurv'});
megakine = innerjoin(maxcurv_tbl, megakine, 'Keys','clip_id');

%% Slip parameter
megakine.slip = megakine.UsBL ./ megakine.wavespeedBLS;

%% Save final kinematics dataset
writetable(megakine, 'Megakine.csv');
writetable(megakine, 'Megakine.xlsx');


%% ---- LOCAL HELPER FUNCTIONS ----

function peakloc = track_peaks(YPeak)
% Tracks wave peaks moving anterior to posterior down the fish body.
% Enforces that peaks only move forward (posterior); if a peak would
% move anteriorly it is held stationary (flagged and removed later).
    clipid   = unique(YPeak.clip_id);
    peakloc  = table();
    
    for i = 1:length(clipid)
        id  = clipid(i);
        dfi = YPeak(YPeak.clip_id == id, :);
        
        frames = unique(dfi.frame);
        min_fr = min(frames);
        
        % Start peak: minimum location at first frame
        pl  = min(dfi.location(dfi.frame == min_fr));
        inc = 1;
        while pl > 150 && inc < height(dfi) % handle clips where first peak is at tail
            inc = inc + 1;
            dfi = dfi(inc:end, :);
            frames = unique(dfi.frame);
            min_fr = min(frames);
            pl = min(dfi.location(dfi.frame == min_fr));
        end
        
        peakloci = table();
        for j = 1:length(frames)
            fr  = frames(j);
            dfj = dfi(dfi.frame == fr, :);
            
            if fr == min_fr
                % Keep initial peak as-is
            else
                % Find next peak posterior to current pl
                candidates = dfj.location(dfj.location > pl);
                if isempty(candidates)
                    % No valid posterior peak: hold position
                else
                    pl = min(candidates);
                end
            end
            peakloci = [peakloci; table(id, fr, pl, ...
                'VariableNames', {'clip_id','frame','location'})];
        end
        peakloc = [peakloc; peakloci];
    end
    
    % Remove rows where peak was held stationary (flag = 1)
    peakloc.flag = zeros(height(peakloc), 1);
    clipids_u = unique(peakloc.clip_id);
    for i = 1:length(clipids_u)
        id   = clipids_u(i);
        mask = find(peakloc.clip_id == id);
        for k = 2:length(mask)
            if peakloc.location(mask(k)) == peakloc.location(mask(k-1))
                peakloc.flag(mask(k)) = 1;
            end
        end
    end
    peakloc = peakloc(peakloc.flag == 0, :);
    peakloc.flag = [];
end


function Midcount = compute_midcount(midlines, Multisp)
% Computes frame range and proportion of tailbeat covered per clip.
    clipids  = unique(midlines.clip_id);
    clip_id  = clipids;
    minframe = nan(length(clipids),1);
    maxframe = nan(length(clipids),1);
    nframes  = nan(length(clipids),1);
    species  = cell(length(clipids),1);
    
    for i = 1:length(clipids)
        id   = clipids(i);
        mask = midlines.clip_id == id;
        minframe(i) = min(midlines.frame(mask));
        maxframe(i) = max(midlines.frame(mask));
        nframes(i)  = maxframe(i) - minframe(i) + 1;
        sp = unique(midlines.species(mask));
        species{i}  = sp{1};
    end
    
    Midcount = table(clip_id, species, minframe, maxframe, nframes);
    proptb   = Multisp(:, {'clip_id','totalframes'});
    Midcount = innerjoin(Midcount, proptb, 'Keys','clip_id');
    Midcount.Prop_TB = Midcount.nframes ./ Midcount.totalframes;
end


function wavemodels = fit_wave_models(peakloc, plot_folder)
% Fits linear regression of PropLoc ~ PropTB for each clip.
% Saves diagnostic plots showing data points and fitted line.
    mkdir(plot_folder);
    clipids    = unique(peakloc.clip_id);
    wavemodels = table();
    
    for i = 1:length(clipids)
        id     = clipids(i);
        a_data = peakloc(peakloc.clip_id == id, :);
        
        coeffs = polyfit(a_data.PropTB, a_data.PropLoc, 1);
        slope  = coeffs(1);
        intcpt = coeffs(2);
        
        y_fit  = polyval(coeffs, a_data.PropTB);
        ss_res = sum((a_data.PropLoc - y_fit).^2);
        ss_tot = sum((a_data.PropLoc - mean(a_data.PropLoc)).^2);
        R2     = 1 - ss_res / ss_tot;
        
        % Standard errors via manual least squares
        n   = height(a_data);
        X   = [ones(n,1), a_data.PropTB];
        b   = X \ a_data.PropLoc;
        res = a_data.PropLoc - X*b;
        s2  = sum(res.^2) / (n-2);
        SE  = sqrt(diag(s2 * inv(X'*X)));
        
        wavemodels = [wavemodels; table(id, intcpt, SE(1), slope, SE(2), R2, ...
            'VariableNames',{'clip_id','Intercept','SE_Intercept','Slope','SE_Slope','R2'})];
        
        % Save plot
        fig = figure('Visible','off');
        plot(a_data.PropTB, a_data.PropLoc, 'ko-'); hold on;
        plot(a_data.PropTB, y_fit, 'b-');
        xlim([0,1]); ylim([-0.02, 1.02]);
        xlabel('Proportion of tailbeat'); ylabel('Proportion of body length');
        title(['Clip ID - ', num2str(id)]);
        saveas(fig, fullfile(plot_folder, [num2str(id), '.png']));
        close(fig);
    end
end