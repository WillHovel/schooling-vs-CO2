%% ANALYSIS - Kinematics + Morphometrics
% Translated from R by E. Goerig

clear; clc;

%% Import datasets
Megakine = readtable('Megakine.csv');
Megakine = Megakine(:, 2:32);

Morpho = readtable('resultsMeansSE.csv');
Morpho = Morpho(:, 2:42);

% Standardize species names (replace spaces with underscores)
Megakine.species = strrep(Megakine.species, ' ', '_');
Morpho.species   = strrep(Morpho.species,   ' ', '_');

% Fix known spelling difference
Morpho.species{strcmp(Morpho.species, 'Danio rerio')} = 'Danio_rerio';

% Check for mismatches between datasets
only_in_megakine = setdiff(unique(Megakine.species), unique(Morpho.species));
only_in_morpho   = setdiff(unique(Morpho.species),   unique(Megakine.species));
disp('In Megakine but not Morpho:'); disp(only_in_megakine);
disp('In Morpho but not Megakine:'); disp(only_in_morpho);

% Merge on species
MegakineMorpho = innerjoin(Megakine, Morpho, 'Keys', 'species');

% Save merged dataset
writetable(MegakineMorpho, 'MegakineMorpho.xlsx');
writetable(MegakineMorpho, 'MegakineMorpho.csv');

%% New Variables

% Fineness ratios - measures of body streamlining
% flateral = SL / max_depth  (lateral view)
% fdorsal  = SL / max_width  (dorsal view)
% f1       = SL / sqrt(max_depth * max_width)  (geometric mean)
% f        = full formula accounting for elliptical cross-section

MegakineMorpho.flateral = MegakineMorpho.standard_length ./ MegakineMorpho.max_depth;
MegakineMorpho.fdorsal  = MegakineMorpho.standard_length ./ MegakineMorpho.max_width;
MegakineMorpho.f1       = MegakineMorpho.standard_length ./ ...
    sqrt(MegakineMorpho.max_depth .* MegakineMorpho.max_width);

rd = sqrt(MegakineMorpho.max_depth);
rb = sqrt(MegakineMorpho.max_width);

% Full fineness ratio formula
MegakineMorpho.f = MegakineMorpho.standard_length ./ ...
    (3*(rd + rb) - sqrt((rd + 3*rb) .* (rb + 3*rd)));

% Specific wavelength (SW) - Nangia et al. 2017
% SW = wavelength / tail amplitude
MegakineMorpho.SW = MegakineMorpho.wavelength ./ MegakineMorpho.tailAmp;

figure; histogram(MegakineMorpho.SW);
title('Specific Wavelength distribution');
disp(summary(MegakineMorpho(:, {'SW'})));


%% Descriptive Statistics

vars_to_summarize = {'headAmp','tailAmp','maxAmp','headtailamp', ...
    'wavelength','maxcurv','UsBL','Re','St'};

for v = 1:length(vars_to_summarize)
    varname = vars_to_summarize{v};
    fprintf('\n--- %s ---\n', varname);
    vals = MegakineMorpho.(varname);
    fprintf('Min: %.4f | Max: %.4f | Mean: %.4f | Median: %.4f | SD: %.4f\n', ...
        min(vals,'omitnan'), max(vals,'omitnan'), ...
        mean(vals,'omitnan'), median(vals,'omitnan'), std(vals,'omitnan'));
    fprintf('IQR: %.4f\n', iqr(vals));
end

% Summary by species and swimming mode
species_list = unique(MegakineMorpho.species);
modes_list   = unique(MegakineMorpho.Swimmingmode);

% Example: headAmp summary by species and swimming mode
headAmp_summary = table();
for s = 1:length(species_list)
    sp   = species_list{s};
    mask = strcmp(MegakineMorpho.species, sp);
    vals = MegakineMorpho.headAmp(mask);
    mode = MegakineMorpho.Swimmingmode(mask);
    mode = mode{1};
    n    = sum(~isnan(vals));
    row  = table({sp},{mode},n, min(vals,'omitnan'), max(vals,'omitnan'), ...
        mean(vals,'omitnan'), median(vals,'omitnan'), std(vals,'omitnan'), ...
        'VariableNames',{'species','Swimmingmode','n','min','max','mean','median','sd'});
    headAmp_summary = [headAmp_summary; row];
end
headAmp_summary = sortrows(headAmp_summary, 'mean');
disp(headAmp_summary);

% 4 model species subset
model_sp = {'Anguilla_rostrata','Scomber_scombrus', ...
            'Salvelinus_fontinalis','Thunnus_albacares'};
mask_model = ismember(MegakineMorpho.species, model_sp);
Modelspecies = MegakineMorpho(mask_model, :);

fprintf('\nModel species wavelength:\n');
disp(summary(Modelspecies(:,{'wavelength'})));

% Save summary tables
writetable(headAmp_summary, 'headamp.xlsx');
writetable(Modelspecies,    'Modelspecies.xlsx');



%% Correlation Matrix

vars = {'maxcurv','UgBL','UsBL','TBF','SL','headAmp','tailAmp','maxAmp', ...
    'headtailamp','Re','St','wavelength','max_width','width90', ...
    'max_depth','depth90','f'};

labels = {'Curvature','Groundspeed','Swimspeed','Tailbeat frequency', ...
    'Stride length','Head amp','Tail amp','Max amp','Head:tail amp', ...
    'Reynolds','Strouhal','Wavelength','Body width','Peduncle width', ...
    'Body depth','Peduncle depth','Fineness ratio'};

% Extract numeric matrix
data_mat = zeros(height(MegakineMorpho), length(vars));
for v = 1:length(vars)
    data_mat(:,v) = MegakineMorpho.(vars{v});
end

% Remove rows with any NaN
valid = all(~isnan(data_mat), 2);
data_mat = data_mat(valid, :);

% Compute correlation matrix and p-values
[R, P] = corrcoef(data_mat);

% Mask non-significant correlations (p >= 0.05)
R_masked = R;
R_masked(P >= 0.05) = 0;

% Plot
figure('Position', [100 100 900 800]);
imagesc(R_masked);
colormap(redblue(256)); % use custom redblue colormap (see helper below)
colorbar;
clim([-1 1]);
set(gca, 'XTick', 1:length(labels), 'XTickLabel', labels, ...
         'YTick', 1:length(labels), 'YTickLabel', labels, ...
         'XTickLabelRotation', 45, 'FontSize', 9);
title('Correlation Matrix (p < 0.05 shown)');
saveas(gcf, 'CorrelationMatrix.pdf');



%% Visual Explorations

% Swimming mode ordering and colors
mode_order  = {'Anguilliform','Sub-carangiform','Carangiform','Thunniform'};
mode_colors = [0.647 0 0.149;   % red    (#A50026) - Anguilliform
               0.192 0.212 0.584; % blue  (#313695) - Sub-carangiform
               0.992 0.682 0.380; % orange (#FDAE61) - Carangiform
               0.455 0.678 0.820];% light blue (#74ADD1) - Thunniform

MegakineMorpho.Swimmingmode = categorical(MegakineMorpho.Swimmingmode, mode_order);

% 1. Swim speed vs head amplitude
figure;
scatter(MegakineMorpho.UsBL, MegakineMorpho.headAmp, 30, 'k', 'filled');
xlim([0 30]); ylim([0 0.4]);
xlabel('Swim speed (BL/s)'); ylabel('Head amplitude (BL)');
title('Swim speed vs Head amplitude');
saveas(gcf, 'USBL_headAmp.pdf');

% 2. Reynolds vs Strouhal, colored by swimming mode
figure; hold on;
for m = 1:length(mode_order)
    mask = MegakineMorpho.Swimmingmode == mode_order{m};
    scatter(MegakineMorpho.Re(mask), MegakineMorpho.St(mask), ...
        40, mode_colors(m,:), 'filled', 'DisplayName', mode_order{m});
end
xlabel('Reynolds number'); ylabel('Strouhal number');
legend('Location','best'); box off;
title('Reynolds vs Strouhal');
saveas(gcf, 'FigureX2_StRe.pdf');

% 3. Effect of swim speed on multiple kinematics variables
kine_vars   = {'TBF','headAmp','maxAmp','wavelength','maxcurv'};
kine_labels = {'Tailbeat frequency (Hz)','Head amplitude (BL)', ...
    'Max amplitude (BL)','Wavelength (BL)','Max curvature (BL)'};

lab_list   = unique(Megakine.lab);
lab_colors = [0 0.545 0.545;   % dark cyan
              0.545 0 0.545;   % dark magenta
              0 0 1;            % blue
              0.855 0.647 0.125]; % dark goldenrod

figure('Position',[100 100 600 900]);
for v = 1:length(kine_vars)
    subplot(length(kine_vars), 1, v); hold on;
    for l = 1:length(lab_list)
        mask = strcmp(Megakine.lab, lab_list{l});
        scatter(Megakine.UsBL(mask), Megakine.(kine_vars{v})(mask), ...
            20, lab_colors(l,:), 'filled', 'DisplayName', lab_list{l});
    end
    ylabel(kine_labels{v}); box off;
    if v == length(kine_vars), xlabel('Swim speed (BL/s)'); end
    if v == 1, legend('Location','best'); end
end
sgtitle('Effect of swim speed on kinematics');
saveas(gcf, 'Swimspeed.pdf');

% 4. Same plot colored by swimming mode
figure('Position',[100 100 600 900]);
for v = 1:length(kine_vars)
    subplot(length(kine_vars), 1, v); hold on;
    for m = 1:length(mode_order)
        mask = MegakineMorpho.Swimmingmode == mode_order{m};
        scatter(MegakineMorpho.UsBL(mask), MegakineMorpho.(kine_vars{v})(mask), ...
            20, mode_colors(m,:), 'filled', 'DisplayName', mode_order{m});
    end
    ylabel(kine_labels{v}); box off;
    if v == length(kine_vars), xlabel('Swim speed (BL/s)'); end
    if v == 1, legend('Location','best'); end
end
sgtitle('Kinematics by swimming mode');
saveas(gcf, 'Figure1_Swimspeed.pdf');



%% PCA Analysis

% --- Helper: run and plot a PCA ---
% (defined as local function at bottom of script)

%% PCA 1 - Full kinematics + morphometrics, averaged by species
vars_pca1 = {'maxcurv','UsBL','TBF','headtailamp','Re','St', ...
    'wavelength','f','max_width','max_depth','depth90','width90'};

[pca1_scores, pca1_coeff, pca1_explained, pca1_labels, pca1_mode] = ...
    run_pca(MegakineMorpho, vars_pca1, mode_order, mode_colors, 'PCA-1');

%% PCA 2 - Ground speed instead of swim speed
vars_pca2 = {'maxcurv','UgBL','TBF','St','maxAmp','headtailamp', ...
    'wavelength','width90','max_depth','max_width','depth90','width90','f'};

[pca2_scores, pca2_coeff, pca2_explained, pca2_labels, pca2_mode] = ...
    run_pca(MegakineMorpho, vars_pca2, mode_order, mode_colors, 'PCA-2');

%% PCA 3 - Reduced set, removes speed/TBF lab effects
vars_pca3 = {'maxcurv','headtailamp','St','wavelength','f', ...
    'max_width','max_depth','depth90','width90'};
labels_pca3 = {'Curvature','Head:tail amp','Strouhal','Wavelength', ...
    'Fineness ratio','Body width','Body depth','Peduncle depth','Peduncle width'};

[pca3_scores, pca3_coeff, pca3_explained, pca3_labels, pca3_mode] = ...
    run_pca(MegakineMorpho, vars_pca3, mode_order, mode_colors, 'PCA-3');

fprintf('PCA3 variance explained by PC1: %.1f%%\n', pca3_explained(1));
fprintf('PCA3 variance explained by PC2: %.1f%%\n', pca3_explained(2));

% Scree plot
figure;
bar(pca3_explained(1:min(9,end)));
xlabel('Principal Component'); ylabel('Variance explained (%)');
title('PCA-3 Scree Plot'); box off;
saveas(gcf, 'Figure4_PCA.pdf');



%% Discriminant Function Analysis (DFA)

dfa_vars = {'maxcurv','maxcurvloc','headAmp','tailAmp','minAmploc','minAmp', ...
    'maxAmploc','maxAmp','f','St','wavelength','depth90','max_depth', ...
    'width90','max_width'};

% Build data matrix
dfa_data = MegakineMorpho(:, [dfa_vars, {'Swimmingmode'}]);
dfa_data = rmmissing(dfa_data);

X = table2array(dfa_data(:, dfa_vars));
y = dfa_data.Swimmingmode;

% Standardize (center + scale, equivalent to R's scale())
X_scaled = (X - mean(X)) ./ std(X);

% Fit LDA using fitcdiscr
lda_model = fitcdiscr(X_scaled, y);

% Cross-validated predictions
cv_model  = crossval(lda_model, 'KFold', 10);
cv_labels = kfoldPredict(cv_model);

% Confusion matrix and accuracy
C = confusionmat(y, cv_labels);
disp('Confusion matrix:'); disp(C);
disp('Per-class accuracy:');
disp(diag(C) ./ sum(C, 2));
fprintf('Overall accuracy: %.1f%%\n', 100 * sum(diag(C)) / sum(C(:)));

% Plot LDA scores (first 2 discriminant dimensions)
[~, scores] = predict(lda_model, X_scaled);
figure; hold on;
for m = 1:length(mode_order)
    mask = strcmp(string(y), mode_order{m});
    scatter(scores(mask,1), scores(mask,2), ...
        40, mode_colors(m,:), 'filled', 'DisplayName', mode_order{m});
end
xlabel('Discriminant Function 1'); ylabel('Discriminant Function 2');
legend('Location','best'); box off;
title('LDA - Swimming mode classification');



%% DBSCAN Density-based Clustering

% Average all variables by species first
species_list = unique(MegakineMorpho.species);
all_vars = {'maxcurv','maxcurvloc','SL','headAmp','tailAmp','minAmp', ...
    'maxAmp','minAmploc','maxAmploc','headtailamp','Re','St','wavelength', ...
    'f','max_width','max_depth','depth90','width90'};

n_sp  = length(species_list);
C_mat = nan(n_sp, length(all_vars));
for s = 1:n_sp
    mask = strcmp(MegakineMorpho.species, species_list{s});
    for v = 1:length(all_vars)
        C_mat(s,v) = mean(MegakineMorpho.(all_vars{v})(mask), 'omitnan');
    end
end

%% DBSCAN on morphometrics
morph_idx = ismember(all_vars, {'f','max_depth','max_width','depth90','width90'});
C_morph   = C_mat(:, morph_idx);
C_morph_s = (C_morph - mean(C_morph)) ./ std(C_morph); % standardize

% Find epsilon using kNN distance plot (equivalent to kNNdistplot)
k = 6;
knn_dist = sort(pdist2(C_morph_s, C_morph_s, 'euclidean'), 2);
knn_k    = sort(knn_dist(:, k+1)); % k+1 because col 1 is self
figure;
plot(knn_k); title('kNN Distance Plot - Morphometrics');
xlabel('Points sorted by distance'); ylabel([num2str(k), '-NN distance']);
% Inspect this plot to choose epsilon visually (elbow point)

epsilon_morph = 0.12; min_pts_morph = 6;
morph_labels  = dbscan(C_morph_s, epsilon_morph, min_pts_morph);

% Pairs plot colored by cluster
cluster_colors = [1 0 0; 0.192 0.212 0.584; 0.992 0.682 0.380]; % red, blue, orange
morph_var_names = {'Fineness ratio','Body depth','Body width','Peduncle depth','Peduncle width'};
plot_pairs(C_morph_s, morph_labels, morph_var_names, cluster_colors, 'B - Morphometrics');
saveas(gcf, 'Figure3A-Morphometrics.pdf');

%% DBSCAN on kinematics
kine_idx  = ismember(all_vars, {'headAmp','tailAmp','wavelength','maxcurv','St'});
C_kine    = C_mat(:, kine_idx);
C_kine_s  = (C_kine - mean(C_kine)) ./ std(C_kine);

knn_dist_k = sort(pdist2(C_kine_s, C_kine_s, 'euclidean'), 2);
knn_k2     = sort(knn_dist_k(:, k+1));
figure;
plot(knn_k2); title('kNN Distance Plot - Kinematics');
xlabel('Points sorted by distance'); ylabel([num2str(k), '-NN distance']);

epsilon_kine = 4.5; min_pts_kine = 6;
kine_labels  = dbscan(C_kine_s, epsilon_kine, min_pts_kine);

kine_var_names = {'Head amplitude','Tail amplitude','Wavelength','Curvature','Strouhal'};
plot_pairs(C_kine_s, kine_labels, kine_var_names, cluster_colors, 'C - Kinematics');
saveas(gcf, 'Figure3B_Kinematics.pdf');



%% ---- LOCAL HELPER FUNCTIONS ----

function [scores, coeff, explained, sp_labels, mode_vec] = ...
        run_pca(data, var_list, mode_order, mode_colors, title_str)
% Averages data by species, runs standardized PCA, and plots biplot.

    species_list = unique(data.species);
    n_sp   = length(species_list);
    X      = nan(n_sp, length(var_list));
    mode_vec = cell(n_sp, 1);

    for s = 1:n_sp
        sp   = species_list{s};
        mask = strcmp(data.species, sp);
        for v = 1:length(var_list)
            X(s,v) = mean(data.(var_list{v})(mask), 'omitnan');
        end
        mode_vec{s} = data.Swimmingmode(find(mask,1));
        if iscell(mode_vec{s}), mode_vec{s} = mode_vec{s}{1}; end
    end

    % Remove species with any NaN
    valid    = all(~isnan(X), 2);
    X        = X(valid, :);
    sp_labels = species_list(valid);
    mode_vec  = mode_vec(valid);

    % Standardize and run PCA
    X_s = (X - mean(X)) ./ std(X);
    [coeff, scores, ~, ~, explained] = pca(X_s);

    fprintf('\n%s - Variance explained:\n', title_str);
    for i = 1:min(4, length(explained))
        fprintf('  PC%d: %.1f%%\n', i, explained(i));
    end

    % Biplot with species labels colored by swimming mode
    figure('Position',[100 100 900 700]); hold on;
    for m = 1:length(mode_order)
        mask_m = strcmp(mode_vec, mode_order{m});
        scatter(scores(mask_m,1), scores(mask_m,2), ...
            60, mode_colors(m,:), 'filled', 'DisplayName', mode_order{m});
        % Label species
        for s = 1:sum(mask_m)
            idx = find(mask_m); 
            text(scores(idx(s),1)+0.1, scores(idx(s),2), ...
                strrep(sp_labels{idx(s)},'_',' '), 'FontSize', 7);
        end
    end

    % Draw loading arrows
    scale = max(abs(scores(:,1:2)),[],'all') * 0.8;
    for v = 1:length(var_list)
        quiver(0, 0, coeff(v,1)*scale, coeff(v,2)*scale, ...
            'k', 'MaxHeadSize', 0.5, 'AutoScale','off');
        text(coeff(v,1)*scale*1.1, coeff(v,2)*scale*1.1, ...
            var_list{v}, 'FontSize', 8, 'Color', [0.3 0.3 0.3]);
    end

    xline(0,'--','Color',[0.7 0.7 0.7]);
    yline(0,'--','Color',[0.7 0.7 0.7]);
    xlabel(sprintf('PC1 (%.1f%% explained)', explained(1)));
    ylabel(sprintf('PC2 (%.1f%% explained)', explained(2)));
    legend('Location','best'); box off;
    title(title_str);
    saveas(gcf, [title_str, '.pdf']);
end


function plot_pairs(X, cluster_labels, var_names, colors, fig_title)
% Creates a pairs/scatterplot matrix colored by DBSCAN cluster.
    n_vars = size(X, 2);
    n_clust = max(cluster_labels) + 1; % +1 because noise = -1 -> 0 after shift
    labels_shifted = cluster_labels + 2; % shift: noise(-1)->1, c1(0)->2, c2(1)->3
    labels_shifted(cluster_labels == -1) = 1; % noise gets color 1

    figure('Position',[100 100 900 800]);
    for r = 1:n_vars
        for c = 1:n_vars
            subplot(n_vars, n_vars, (r-1)*n_vars + c);
            if r == c
                histogram(X(:,r), 10, 'FaceColor',[0.7 0.7 0.7]);
            else
                hold on;
                for cl = 1:size(colors,1)
                    mask = labels_shifted == cl;
                    scatter(X(mask,c), X(mask,r), 20, colors(cl,:), 'filled');
                end
            end
            if r == n_vars, xlabel(var_names{c},'FontSize',7); end
            if c == 1,      ylabel(var_names{r},'FontSize',7); end
            box off; set(gca,'FontSize',6);
        end
    end
    sgtitle(fig_title);
end


function cmap = redblue(n)
% Generates a red-white-blue colormap (equivalent to RdBu in R/corrplot).
    if nargin < 1, n = 256; end
    half = floor(n/2);
    r1 = linspace(0.647, 1, half);
    g1 = linspace(0, 1, half);
    b1 = linspace(0.149, 1, half);
    r2 = linspace(1, 0.192, n-half);
    g2 = linspace(1, 0.212, n-half);
    b2 = linspace(1, 0.584, n-half);
    cmap = [[r1,r2]', [g1,g2]', [b1,b2]'];
end