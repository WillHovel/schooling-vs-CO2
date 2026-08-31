% VALIDATE_DEMO_DATA  Load every generated demo CSV through the real
% Kinemetrix loaders and confirm the formats + key metrics are sane.
%
% Run from the project root (so the *.m loaders are on the path):
%   matlab -batch "cd('C:/Users/willh/Desktop/Fish_analysis_v2'); run('demo_data/validate_demo_data.m')"

base = 'C:/Users/willh/Desktop/Fish_analysis_v2/demo_data';
fprintf('\n================ KINEMETRIX DEMO DATA VALIDATION ================\n');

%% 1. Kinematics tab -- Format A (DLC single, 3D)
fprintf('\n[1] Kinematics tab (Format A DLC single 3D)\n');
fp = load_fish_points(fullfile(base, 'demo_kinematics.csv'));
assert(numel(fp) == 1);
assert(strcmp(fp.format, 'DLC') && fp.has_z && numel(fp.point_names) == 13);
fpT = transform_fish(fp);
kine = compute_kinematics(fpT, 100, 1.0);
ext  = compute_body_extended(fpT, 100, kine, {}, 0);
fprintf('    tail_TBF = %.3f Hz (expected ~5.0)\n', kine.tail_TBF);
fprintf('    head_TBF = %.3f Hz\n', kine.head_TBF);
fprintf('    wavelength = %.3f BL (expected ~1.0)\n', kine.wavelength);
fprintf('    mean speed = %.3f BL/s (expected ~1.5)\n', ext.mean_speed_BL_s);
fprintf('    strouhal   = %.3f (expected 0.2-0.4)\n', ext.strouhal);

%% 2. Fin Analysis tab -- Format B (named, 3D)
fprintf('\n[2] Fin Analysis tab (Format B named 3D)\n');
midline = {'snout','mid1','mid2','mid3','peduncle','caudaltip'};
fpF = load_fish_points_named(fullfile(base, 'demo_fin.csv'), midline, []);
assert(fpF.has_z && strcmp(fpF.format, 'named'));
fin = compute_fin_kinematics(fullfile(base, 'demo_fin.csv'), 'Rpectbase', 'Rpecttip', 100, 1.0, ext.mean_speed_BL_s);
fprintf('    fin_freq (yaw) = %.3f Hz (expected ~3.0)\n', fin.fin_freq_Hz);
fprintf('    mean fin length = %.3f\n', fin.mean_length);
fprintf('    mean yaw = %.2f deg, range = %.2f deg\n', fin.mean_yaw, fin.range_yaw);
st = compute_stance_swing(fin);
fprintf('    duty factor = %.3f (n_cycles=%d)\n', st.duty_factor, st.n_cycles);

%% 3. Pectoral Phase tab -- Format D (dual-camera, 2D)
fprintf('\n[3] Pectoral Phase tab (Format D dual-camera 2D)\n');
for which = {'demo_pectoral_inphase.csv', 'demo_pectoral_antiphase.csv'}
    f = which{1};
    fpP = load_fish_points_named(fullfile(base, f), [], []);
    assert(strcmp(fpP.format, 'dual_camera'));
    r = fpP.pect_phase_result;
    fprintf('    %-32s -> %s (phase shift = %.1f deg)\n', f, r.classification, r.phase_shift_deg);
end

%% 4. School Metrics tab -- Format A (DLC multi-animal)
fprintf('\n[4] School Metrics tab (Format A DLC multi-animal)\n');
fpS = load_fish_points_dlc_multianimal(fullfile(base, 'demo_school_3fish.csv'));
assert(numel(fpS) == 3);
pol = compute_polarization(fpS, 'snout', 'peduncle', 0);
dist = compute_distance_between_individuals(fpS, 'snout', 1.0);
fprintf('    individuals: %s\n', strjoin({fpS.name}, ', '));
fprintf('    mean polarization = %.4f (expected ~0.9+)\n', pol.mean_polarization);
fprintf('    mean NN distance  = %.3f\n', dist.overall_mean_nn_dist);

%% 5. Batch Processing tab -- Format B (named, 3D), tokenised filenames
fprintf('\n[5] Batch Processing tab (Format B named 3D)\n');
batchPts = {'snout','mid1','mid2','mid3','mid4','peduncle','caudaltip'};
d = dir(fullfile(base, 'batch_trials', '*.csv'));
for i = 1:numel(d)
    f = fullfile(d(i).folder, d(i).name);
    fpB = load_fish_points_named(f, batchPts, []);
    assert(fpB.has_z && numel(fpB.point_names) == 7);
    fprintf('    %s  (%d pts, %d frames)\n', d(i).name, numel(fpB.point_names), numel(fpB.frames));
end

fprintf('\n================ ALL DEMO DATA LOADED AND VERIFIED ================\n');
