function test_compute_kinematics()
% TEST_COMPUTE_KINEMATICS  Ground-truth checks for compute_kinematics.m.
%
%   All ground truth comes from synth_fish (exact known camera geometry):
%    1. straight horizontal fish, lambda=1, constant A: the fully exact
%       case: TBFs land exactly on the FFT bin, Y_interp must reproduce
%       the known wave (after the documented spatial-mean removal) to 1e-9,
%       amplitudes equal the known time-mean (2/pi)*A/L, wavelength = 1,
%       maxCurv = the known wave curvature, wave_speed = wavelength*tail_TBF;
%    2. lambda = 0.8 (non-integer, 201 raw stations -> spline branch);
%    3. modulated amplitude envelope A(s): the full amplitude profile,
%       head/tail ratio, and max/min locations match the analytic envelope;
%    4. tilted path (30 deg): the case that carries the line-wobble
%       contamination; if the pipeline is exact, wavelength still
%       measures 1 to within 2%;
%    5. circling fish: TBFs and wavelength survive a turning path;
%    6. standing wave: the wavelength gate must return NaN, not a
%       fabricated number;
%    7. min_freq above the signal: everything NaN;
%    8. missing frames: NaN rows stay NaN, validity count exact, metrics
%       survive a 6-frame gap;
%    9. 3-D pitch: Z amplitudes and Z beat frequencies (incl. headZ_TBF,
%       measured at the station behind the head);
%   10. all-NaN fish: zero usable frames, all-NaN outputs, no crash.

    fps = 100; nF = 400; f0 = 2; A = 0.05; L = 10; A_bl = A/L;
    min_freq = 0.5;
    amp_expected = (2/pi) * A_bl;              % time-mean |wave| = 0.003183
    maxCurv_expected = (2/pi) * A_bl * (2*pi)^2;   % time-mean max |curvature| = 0.1257

    % ---- 1) Straight horizontal, lambda=1: the fully exact case ----
    fp = synth_fish('fps', fps, 'nFrames', nF, 'f0', f0, 'lambda', 1, ...
                    'A', A, 'U', 1, 'L', L, 'heading', 0);
    tf = transform_fish(fp);
    k = compute_kinematics(tf, fps, min_freq);

    check(k.n_frames_with_data == nF, 'all frames usable', 'got %d', k.n_frames_with_data);
    check_close(k.head_TBF, f0, 1e-9, 'head_TBF = f0');
    check_close(k.tail_TBF, f0, 1e-9, 'tail_TBF = f0');
    check_close(k.head_fft_freq, f0, 1e-9, 'head_fft_freq = f0');
    check(any(k.head_fft_power > 0), 'head FFT power > 0');
    check_close(k.spline_freq_Hz, f0, 1e-9, 'spline_freq = f0');

    % Y_interp must be the exact known wave minus its per-frame spatial mean
    t = fp.gt.t; s = k.s_norm;
    wave = A_bl * sin(2*pi*(f0*t - s));
    wave = wave - mean(wave, 2);
    check_all_close(k.Y_interp, wave, 1e-9, 'Y_interp == known wave (mean-removed)');
    check_all_close(k.X_interp, repmat(s, nF, 1), 1e-9, 'X_interp = station grid');

    check_close(k.headAmp, amp_expected, 2e-2, 'headAmp = (2/pi) A/L');
    check_close(k.tailAmp, amp_expected, 2e-2, 'tailAmp = (2/pi) A/L');
    check_close(k.headTailAmpRatio, 1, 2e-2, 'constant A: head/tail ratio = 1');
    check_close(k.minAmp, amp_expected, 2e-2, 'minAmp = uniform amplitude');
    check_close(k.maxAmp, amp_expected, 2e-2, 'maxAmp = uniform amplitude');

    check_close(k.wavelength, 1, 1e-2, 'wavelength = 1');
    check_close(k.wave_spatial_freq, 1, 1e-2, 'spatial freq = 1/BL');
    check(k.wave_power >= 0.99, 'wavelength fit r2 >= 0.99', 'got %.4f', k.wave_power);
    check_close(k.wave_speed_BL_s, f0, 1e-2, 'wave_speed = 2 BL/s');
    check_all_close(k.wave_speed_BL_s, k.wavelength * k.tail_TBF, 1e-12, ...
                    'wave_speed identity');

    check_close(k.maxCurv, maxCurv_expected, 2e-2, 'maxCurv = (2/pi) A/L (2pi)^2');
    check(isempty(k.curv3d_mean), 'no 3-D curvature in 2-D data');

    % ---- 2) lambda = 0.8, non-integer (spline interpolation branch) ----
    fp8 = synth_fish('fps', fps, 'nFrames', nF, 'nPts', 201, 'f0', f0, 'lambda', 0.8, ...
                     'A', A, 'U', 1, 'L', L, 'heading', 0);
    tf8 = transform_fish(fp8);
    k8 = compute_kinematics(tf8, fps, min_freq);
    check_close(k8.head_TBF, f0, 1e-9, 'lambda=0.8: head_TBF');
    check_close(k8.tail_TBF, f0, 1e-9, 'lambda=0.8: tail_TBF');
    check_close(k8.wavelength, 0.8, 1e-2, 'lambda=0.8: wavelength measured');
    check(k8.wave_power >= 0.99, 'lambda=0.8: fit r2', 'got %.4f', k8.wave_power);
    check_close(k8.wave_speed_BL_s, 0.8*f0, 1e-2, 'lambda=0.8: wave_speed');

    % ---- 3) Modulated amplitude envelope A(s) = A0*(1 + 0.5 sin(2 pi s)) ----
    Am = @(ss) A .* (1 + 0.5 .* sin(2*pi*ss));
    fpm = synth_fish('fps', fps, 'nFrames', nF, 'f0', f0, 'lambda', 1, ...
                     'A', Am, 'U', 1, 'L', L, 'heading', 0);
    tfm = transform_fish(fpm);
    km = compute_kinematics(tfm, fps, min_freq);

    sm = km.s_norm;
    G = (Am(sm)/L) .* exp(-2i*pi*sm);
    M = mean(G);
    profile = (2/pi) * abs(G - M);
    [p_max, i_max] = max(profile);
    [p_min, i_min] = min(profile);
    head_sel = sm <= 0.05; tail_sel = sm >= 0.95;
    headAmp_p = mean(profile(head_sel));
    tailAmp_p = mean(profile(tail_sel));

    check_all_close(km.amp_mean, profile, 2.5e-4, 'modulated amplitude profile');
    check_close(km.maxAmp, p_max, 4e-2, 'modulated maxAmp');
    check_close(km.minAmp, p_min, 4e-2, 'modulated minAmp');
    check_close(km.maxAmpLoc, sm(i_max), 2e-2, 'modulated maxAmpLoc (~0.25)');
    check_close(km.minAmpLoc, sm(i_min), 2e-2, 'modulated minAmpLoc (~0.75)');
    check_close(km.headAmp, headAmp_p, 4e-2, 'modulated headAmp');
    check_close(km.tailAmp, tailAmp_p, 4e-2, 'modulated tailAmp');
    check_close(km.headTailAmpRatio, headAmp_p/tailAmp_p, 6e-2, 'modulated head/tail ratio');
    check_close(km.wavelength, 1, 1e-2, 'modulated: wavelength unaffected by envelope');
    % At the envelope peak (s = 0.25, E = 1.5 A) the second derivative is
    % (E'' - E*(2 pi)^2)*sin(phi) with E'' = -0.5 A (2 pi)^2, so the
    % curvature coefficient is 2*A*(2 pi)^2, not 1.5*A*(2 pi)^2 (the
    % envelope's own curvature adds to the wave's).
    check_close(km.maxCurv, (2/pi) * 2 * A_bl * (2*pi)^2, 4e-2, 'modulated maxCurv');
    check_close(km.maxCurvLoc, 0.25, 3e-2, 'modulated maxCurvLoc = amplitude peak');

    % ---- 4) Tilted path (30 deg): line-wobble contamination gate ----
    fpT = synth_fish('fps', fps, 'nFrames', nF, 'nPts', 201, 'f0', f0, 'lambda', 1, ...
                     'A', A, 'U', 1, 'L', L, 'heading', deg2rad(30));
    tfT = transform_fish(fpT);
    kT = compute_kinematics(tfT, fps, min_freq);
    check_close(kT.head_TBF, f0, 5e-2, 'tilted head_TBF');
    check_close(kT.tail_TBF, f0, 5e-2, 'tilted tail_TBF');
    fprintf('  [tilted-path wavelength measured = %.4f BL (true 1.0)]\n', kT.wavelength);
    check_close(kT.wavelength, 1, 2e-2, 'tilted-path wavelength (line-wobble gate)');
    check(kT.wave_power >= 0.9, 'tilted fit r2', 'got %.4f', kT.wave_power);

    % ---- 5) Circling fish ----
    % Turn rate 0.01 rad/s -> circle radius U/omega = 100 (10 BL). A fast
    % turn (omega = U/R with R ~ 0.5 BL) is degenerate: the body-frame
    % lateral swing and the rotating-heading s-linear projection dwarf the
    % wave, which the constant-heading contamination model cannot absorb.
    % 201 raw stations like the tilted case: the s-linear camera-projection
    % ramp in the corrected signal rings under fft_interp's band limit, and
    % at the 21-point default (10 cycles/BL) the Gibbs ripple swamps the
    % curvature metric.
    fpc = synth_fish('fps', fps, 'nFrames', nF, 'nPts', 201, 'f0', f0, 'lambda', 1, 'A', A, ...
                     'U', 1, 'L', L, 'path', 'circle', 'heading', @(tt) 0.01*tt);
    tfc = transform_fish(fpc);
    kc = compute_kinematics(tfc, fps, min_freq);
    check_close(kc.head_TBF, f0, 5e-2, 'circling head_TBF');
    check_close(kc.tail_TBF, f0, 5e-2, 'circling tail_TBF');
    check_close(kc.wavelength, 1, 1e-1, 'circling wavelength');
    check(kc.wave_power >= 0.6, 'circling fit r2', 'got %.4f', kc.wave_power);
    check_close(kc.maxCurv, maxCurv_expected, 1.5e-1, 'circling maxCurv');

    % ---- 6) Standing wave: wavelength gate must produce NaN ----
    % lambda = inf makes the wave EXACTLY standing: every station shares
    % one oscillation, so the spatial-mean removal zeroes the signal to
    % machine precision and the no-oscillation gate must fire. (A finite
    % lambda like 1e6 leaves an s-varying residual of A*2*pi/lambda that
    % survives the gate and fits as a garbage wavelength.)
    fpsw = synth_fish('fps', fps, 'nFrames', nF, 'f0', f0, 'lambda', inf, ...
                      'A', A, 'U', 1, 'L', L, 'heading', 0);
    tfsw = transform_fish(fpsw);
    ksw = compute_kinematics(tfsw, fps, min_freq);
    check_close(ksw.head_TBF, f0, 5e-2, 'standing head_TBF');
    check(isnan(ksw.wavelength), 'standing wave: wavelength NaN (flat-gradient gate)');
    check(isempty(ksw.wave_power), 'standing wave: no fit stats exported');
    check(isnan(ksw.wave_speed_BL_s), 'standing wave: wave_speed NaN');
    check_close(ksw.headAmp, 0, 1e-6, 'standing wave: mean removal leaves no lateral signal');

    % ---- 7) min_freq above the signal: everything NaN ----
    kg = compute_kinematics(tf, fps, 3);
    check(isnan(kg.head_TBF) && isnan(kg.tail_TBF) && isnan(kg.spline_freq_Hz), ...
          'min_freq gate: TBFs NaN');
    check(isnan(kg.wavelength) && isnan(kg.wave_speed_BL_s), 'min_freq gate: wavelength NaN');

    % ---- 8) Missing frames ----
    fpmiss = fp;
    fpmiss.points(10:12, 1, :) = NaN;      % head missing on 3 frames
    fpmiss.points(20:22, end, :) = NaN;    % tail missing on 3 frames
    tfmiss = transform_fish(fpmiss);
    kmiss = compute_kinematics(tfmiss, fps, min_freq);
    check(kmiss.n_frames_with_data == nF - 6, 'missing frames counted', ...
          'got %d', kmiss.n_frames_with_data);
    check(all(isnan(kmiss.Y_interp(10:12, :)), 'all'), 'missing-frame rows stay NaN');
    check_close(kmiss.head_TBF, f0, 5e-2, 'TBF survives 6-frame gap');
    check_close(kmiss.wavelength, 1, 1e-2, 'wavelength survives 6-frame gap');

    % ---- 9) 3-D pitch oscillation ----
    nPts3 = 201;
    s3 = linspace(0, 1, nPts3);
    zA = 0.02;
    fpz.name = 'pitch_synth';
    fpz.X = repmat(s3, nF, 1);
    fpz.Y = A_bl * sin(2*pi*(f0*t - s3));
    fpz.Z = (zA/L) * sin(2*pi*f0*t) * s3;
    fpz.has_z = true;
    fpz.pre_transformed = true;
    kz = compute_kinematics(fpz, fps, min_freq);
    check_close(kz.headZ_TBF, f0, 1e-9, 'headZ_TBF = f0 (station behind head)');
    check_close(kz.tailZ_TBF, f0, 1e-9, 'tailZ_TBF = f0');
    zAmp = (2/pi) * zA / L;
    check_close(kz.maxAmpZ, zAmp, 2e-2, 'maxAmpZ = (2/pi) zA/L');
    check_close(kz.maxAmpZLoc, 1, 1e-9, 'maxAmpZLoc = tail');
    check_close(kz.minAmpZ, 0, 1e-6, 'minAmpZ = 0 at head');
    check_close(kz.minAmpZLoc, 0, 1e-9, 'minAmpZLoc = head');
    check_close(kz.headAmpZ, zAmp * mean(kz.s_norm(head_sel)), 2e-2, 'headAmpZ');
    check_close(kz.tailAmpZ, zAmp * mean(kz.s_norm(tail_sel)), 2e-2, 'tailAmpZ');
    check_close(kz.maxCurv3D, maxCurv_expected, 5e-2, '3-D curvature ~ 2-D');

    % ---- 10) All-NaN fish: no crash, all-NaN outputs ----
    fpz2 = synth_fish('fps', fps, 'nFrames', nF, 'A', A, 'U', 1, 'L', L, 'heading', 0);
    fpz2.points(:) = NaN;
    tfz = transform_fish(fpz2);
    kz2 = compute_kinematics(tfz, fps, min_freq);
    check(kz2.n_frames_with_data == 0, 'zero-data count', 'got %d', kz2.n_frames_with_data);
    check(isnan(kz2.head_TBF) && isnan(kz2.wavelength) && isnan(kz2.maxCurv), ...
          'zero-data outputs all NaN');
end
