function test_fin_kinematics()
% TEST_FIN_KINEMATICS  Ground-truth checks for compute_fin_kinematics.m
%                      and compute_stance_swing.m.
%
%   Fin: root pinned at the origin, tip on the unit circle at a known
%   angle th(t) = 30 deg * sin(2*pi*2*t). Then yaw == th exactly,
%   fin_length == 1 exactly, fin_freq = 2 Hz (exact FFT bin), stride =
%   body_speed/freq, total_dist = the known chord-sum of the tip's arc
%   (2*|sin(dTh/2)| summed, the exact identity for unit-circle steps).
%   Compound roots ('A+B'), all-NA tip columns, and missing columns are
%   exercised too.
%
%   Stance/swing: a hand-built fin struct with 10 cycles of 30 slow
%   (stance) + 15 fast (swing) frames. p95 = 100, threshold = 25, so the
%   classification is exact. Frame 1 has no prior frame, so the first
%   stance bout is [2 30]: mean contact = 0.299 s, swing = 0.15 s, duty =
%   0.299/0.449. A 1-frame tracking gap inside a stance bout must be
%   bridged (2 unknown frames <= MaxGapFrames) without splitting the bout.

    d = tempname;
    mkdir(d);
    cleanup = onCleanup(@() rmdir(d, 's'));
    f = @(nm) fullfile(d, nm);

    fps = 100; nF = 100;
    t = (0:nF-1)'/fps;
    th = deg2rad(30) * sin(2*pi*2*t);          % 2 Hz yaw oscillation, no wrap
    root = zeros(nF, 3);
    tip  = [cos(th), sin(th), zeros(nF,1)];
    rows = [ {'root_X','root_Y','root_Z','tip_X','tip_Y','tip_Z'}
             num2cell([root, tip]) ];
    write_csv(f('fin.csv'), rows);

    fin = compute_fin_kinematics(f('fin.csv'), 'root', 'tip', fps, 0.5, 0.1);
    check_close(fin.fin_freq_Hz, 2, 1e-9, 'fin_freq = 2 Hz');
    check_close(fin.mean_length, 1, 1e-9, 'fin length = 1');
    check_close(fin.std_length, 0, 1e-9, 'length constant');
    check_all_close(fin.yaw, rad2deg(th), 1e-7, 'yaw = known angle series');
    check_close(fin.mean_yaw, 0, 1e-3, 'mean yaw = 0');
    expected_range = 60 * max(sin(2*pi*2*t));    % sampled grid misses the exact peak
    check_close(fin.range_yaw, expected_range, 1e-7, 'yaw range = 60 deg on the sampled grid');
    check(isnan(fin.fin_freq_pitch_Hz), 'no pitch motion -> pitch freq NaN');
    check_close(fin.stride_duration_s, 0.5, 1e-9, 'stride duration = 0.5 s');
    check_close(fin.stride_length_BL, 0.05, 1e-9, 'stride = body_speed/freq');
    check(fin.n_valid == nF, 'all frames valid', 'got %d', fin.n_valid);
    check(fin.pct_valid == 100, 'pct_valid = 100', 'got %.1f', fin.pct_valid);

    expected_total = sum(2*abs(sin(diff(th)/2)));   % exact unit-circle arc
    check_close(fin.total_dist, expected_total, 1e-8, 'total_dist = known arc');
    check_close(fin.mean_speed, expected_total/(nF-1)*fps, 1e-8, 'mean tip speed');

    % ---- Compound root = mean of pair ----
    rows2 = [ {'root1_X','root1_Y','root1_Z','root2_X','root2_Y','root2_Z', ...
               'tip_X','tip_Y','tip_Z'}
              num2cell([root, root, tip]) ];
    write_csv(f('fin_comp.csv'), rows2);
    finC = compute_fin_kinematics(f('fin_comp.csv'), 'root1+root2', 'tip', fps, 0.5);
    check_all_close(finC.root_xyz, root, 1e-12, 'compound root = mean of pair');
    check_close(finC.fin_freq_Hz, 2, 1e-9, 'compound-root freq = 2 Hz');

    % ---- All-NA tip column -> clean error, not a crash ----
    rows3 = rows;
    for r = 2:size(rows3,1)
        rows3{r,4} = 'NA'; rows3{r,5} = 'NA'; rows3{r,6} = 'NA';
    end
    write_csv(f('fin_na.csv'), rows3);
    threw = false;
    try
        compute_fin_kinematics(f('fin_na.csv'), 'root', 'tip', fps);
    catch err
        threw = true;
        check(contains(err.message, 'fewer than 2 valid frames'), ...
              'all-NA tip errors with valid-frame message', 'msg = %s', err.message);
    end
    check(threw, 'all-NA tip must error');

    % ---- Missing column -> error naming the column ----
    threw = false;
    try
        compute_fin_kinematics(f('fin.csv'), 'root', 'missing', fps);
    catch err
        threw = true;
        check(contains(err.message, 'not found'), 'missing column error', ...
              'msg = %s', err.message);
    end
    check(threw, 'missing column must error');

    % ==================== STANCE / SWING ====================
    nSt = 450;
    speed = zeros(nSt, 1);
    for c = 0:9
        i0 = c*45 + 1;
        speed(i0:i0+29) = 1;      % stance: 30 slow frames
        speed(i0+30:i0+44) = 100; % swing: 15 fast frames
    end
    finS = struct('tip_speed', speed, 'fps', fps, 'valid', true(nSt,1), ...
                  'rootName', 'r', 'tipName', 't');
    st = compute_stance_swing(finS);
    check_close(st.threshold_used, 25, 1e-9, 'threshold = 0.25*p95');
    check_close(st.mean_contact_time_s, 0.299, 1e-9, 'contact = 0.299 s');
    check_close(st.mean_swing_time_s, 0.15, 1e-9, 'swing = 0.15 s');
    check_close(st.duty_factor, 0.299/0.449, 1e-9, 'duty factor');
    check(st.n_cycles == 10, '10 cycles found', 'got %d', st.n_cycles);
    check(isequal(st.stance_bouts_frames(1,:), [2 30]), 'first stance bout [2 30]');
    check(isequal(st.swing_bouts_frames(1,:), [31 45]), 'first swing bout [31 45]');
    check(st.is_stance(20) && ~st.is_stance(40), 'phase classification exact');
    check_close(st.pct_valid_consecutive, 100*449/450, 1e-9, 'pct consecutive known');

    % ---- 1-frame tracking gap inside stance must be bridged ----
    finS2 = finS; finS2.valid(50) = false;
    st2 = compute_stance_swing(finS2);
    check(size(st2.stance_bouts_frames,1) == 10, 'gap bridged: still 10 stance bouts', ...
          'got %d', size(st2.stance_bouts_frames,1));
    check(isequal(st2.stance_bouts_frames(2,:), [46 75]), ...
          '2 unknown frames inside bout bridged, not split');
end
