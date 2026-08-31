function test_body_extended()
% TEST_BODY_EXTENDED  Ground-truth checks for compute_body_extended.m.
%
%   1. Straight (A=0) 3D fish, heading 30 deg, z_head = 0.02*t, linear
%      z_off ramp: body_angle must be exactly 30 deg with zero std/range/
%      angular velocity, speed exactly U/L = 0.1 BL/s, head_pitch exactly
%      atand((z1-z2)/horiz) = atand(-0.05), head_Z_raw = 0.02*t.
%   2. Undulating fish (A=0.05, lambda=1, f0=2): mean speed still
%      0.1 BL/s, stride length = speed/tail_TBF = 0.05, tail_amp_pp =
%      2*A/L = 0.01, Strouhal = TBF*A_pp/U = 0.2.
%   3. Signed flow correction: flow +0.05 -> through-water speed 0.15
%      ('against'); flow -0.02 -> 0.08 ('with'), exact on a straight fish.
%   4. Hand-built roll fish: L/R pectoral pair tilted 20 deg about the
%      body axis -> roll = 20 deg exactly, level head -> pitch 0.
%   5. Pre-transformed (CURVES) fish translating at 0.1 BL/s: speed exact
%      (bl = 1 by definition), body angle stays NaN.

    fps = 100; nF = 400;

    % ---- 1) Straight 3D fish, heading 30 deg, known z motion ----
    fp = synth_fish('fps', fps, 'nFrames', nF, 'A', 0, 'U', 1, 'L', 10, ...
                    'heading', deg2rad(30), 'z_head', @(tt) 0.02*tt, ...
                    'z_off', 0.5*linspace(0, 1, 21));
    tf = transform_fish(fp);
    ext = compute_body_extended(tf, fps, [], [], 0);
    t = fp.gt.t;

    check_all_close(ext.body_angle_deg, 30*ones(nF,1), 1e-9, 'body_angle = 30 deg');
    check_close(ext.mean_body_angle_deg, 30, 1e-9, 'mean body angle = 30');
    check_close(ext.std_body_angle_deg, 0, 1e-9, 'std body angle = 0');
    check_close(ext.range_body_angle_deg, 0, 1e-9, 'range body angle = 0');
    check(isnan(ext.angular_velocity_deg_s(1)), 'first-frame angular velocity NaN');
    check_all_close(ext.angular_velocity_deg_s(2:end), zeros(nF-1,1), 1e-9, 'no turning');
    check_close(ext.mean_speed_BL_s, 0.1, 1e-9, 'speed = U/L = 0.1 BL/s');
    check_close(ext.std_speed_BL_s, 0, 1e-9, 'std speed = 0');
    check_close(ext.peak_speed_BL_s, 0.1, 1e-9, 'peak speed = 0.1');
    check_all_close(ext.head_pitch_deg, atand(-0.05)*ones(nF,1), 1e-9, ...
                    'head_pitch = atand(-0.025/0.5)');
    check_close(ext.mean_head_pitch_deg, atand(-0.05), 1e-9, 'mean head pitch');
    check_close(ext.range_head_pitch_deg, 0, 1e-9, 'pitch range = 0');
    check_all_close(ext.head_Z_raw, 0.02*t, 1e-9, 'head_Z_raw = 0.02 t');
    check(isnan(ext.stride_length_BL), 'no kine -> stride NaN');
    check(isnan(ext.strouhal), 'no kine -> strouhal NaN');
    check(strcmp(ext.flow_orientation, 'none'), 'flow orientation "none"');
    check_close(ext.tail_amp_pp_BL, 0, 1e-12, 'straight fish tail amp = 0');
    check(~ext.roll_available, 'no roll pair -> roll unavailable');
    check(all(isnan(ext.roll_deg)), 'roll NaN without pair');

    % ---- 2) Undulating fish: stride / tail amp / Strouhal ----
    fp2 = synth_fish('fps', fps, 'nFrames', nF, 'f0', 2, 'lambda', 1, 'A', 0.05, ...
                     'U', 1, 'L', 10, 'heading', 0);
    tf2 = transform_fish(fp2);
    kine2 = compute_kinematics(tf2, fps, 0.5);
    ext2 = compute_body_extended(tf2, fps, kine2, [], 0);
    check_close(ext2.mean_speed_BL_s, 0.1, 1e-3, 'undulating mean speed');
    check_close(ext2.stride_length_BL, 0.05, 2e-2, 'stride = speed/tail_TBF');
    check_close(ext2.tail_amp_pp_BL, 0.01, 1e-1, 'tail amp pp = 2*A/L');
    check_close(ext2.strouhal, 0.2, 1e-1, 'St = TBF*A_pp/U');

    % ---- 3) Signed flow correction (exact on the straight fish) ----
    extA = compute_body_extended(tf, fps, [], [], 0.05);
    check(strcmp(extA.flow_orientation, 'against'), 'flow +0.05 -> against');
    check_close(extA.flow_BL_s, 0.05, 0, 'flow stored');
    check_close(extA.mean_speed_through_water_BL_s, 0.15, 1e-9, 'U_tw = 0.1 + 0.05');
    extW = compute_body_extended(tf, fps, [], [], -0.02);
    check(strcmp(extW.flow_orientation, 'with'), 'flow -0.02 -> with');
    check_close(extW.mean_speed_through_water_BL_s, 0.08, 1e-9, 'U_tw = |0.1 - 0.02|');

    % ---- 4) Roll: hand-built fish with tilted L/R pectoral pair ----
    nFr = 20; t20 = tand(20);
    pts4 = zeros(nFr, 5, 3);
    pts4(:, 1, :) = repmat([0 0 0],     nFr, 1);   % snout
    pts4(:, 2, :) = repmat([0.5 0 0],   nFr, 1);   % mid
    pts4(:, 3, :) = repmat([1 0 0],     nFr, 1);   % tail
    pts4(:, 4, :) = repmat([0.5 -0.5 0],    nFr, 1);   % LPect
    pts4(:, 5, :) = repmat([0.5 0.5 t20],   nFr, 1);   % RPect
    fpr.name = 'roll_synth';
    fpr.frames = (1:nFr)';
    fpr.point_names = {'snout','mid','tail','LPect','RPect'};
    fpr.points = pts4;
    fpr.has_z = true;
    fpr.bl_per_frame = ones(nFr, 1);
    extR = compute_body_extended(fpr, 100, [], {'LPect','RPect'}, 0);
    check(extR.roll_available, 'roll available with pair');
    check_all_close(extR.roll_deg, 20*ones(nFr,1), 1e-9, 'roll = atan2d(tan20, 1) = 20 deg');
    check_close(extR.mean_roll_deg, 20, 1e-9, 'mean roll = 20');
    check_close(extR.std_roll_deg, 0, 1e-9, 'std roll = 0');
    check_all_close(extR.head_pitch_deg, zeros(nFr,1), 1e-9, 'level head -> pitch 0');
    check_close(extR.mean_speed_BL_s, 0, 1e-9, 'static fish speed = 0');

    % ---- 5) Pre-transformed (CURVES): speed in BL, angle NaN ----
    nF5 = 100;
    x5 = linspace(0, 1, 11);
    pts5 = zeros(nF5, 11, 2);
    pts5(:, :, 1) = x5 + 0.001*(0:nF5-1)';
    fpe.name = 'curves_synth';
    fpe.frames = (1:nF5)';
    fpe.point_names = arrayfun(@(k) sprintf('s%02d', k), 1:11, 'UniformOutput', false);
    fpe.points = pts5;
    fpe.has_z = false;
    fpe.pre_transformed = true;
    extP = compute_body_extended(fpe, 100, [], [], 0);
    check_close(extP.mean_speed_BL_s, 0.1, 1e-9, 'pre-transformed speed = 0.1 BL/s');
    check(all(isnan(extP.body_angle_deg)), 'pre-transformed body angle NaN');
    check(isnan(extP.mean_body_angle_deg), 'pre-transformed mean angle NaN');
end
