function test_filter_dlc_jumps()
% TEST_FILTER_DLC_JUMPS  Ground-truth checks for filter_dlc_jumps.m.
%
%   A straight (A=0), constant-velocity fish has the useful property that
%   the local median of any symmetric neighbor window is EXACTLY the true
%   position at the tested frame (median of evenly spaced velocity offsets
%   is 0), so:
%     1. clean data must produce ZERO filtered point-frames;
%     2. a single-frame teleport must be caught on exactly that frame;
%     3. a 5-consecutive-frame teleport with the default window (13) must
%        catch exactly those 5 frames with no false positives on neighbors
%        (the documented validated behavior);
%     4. a jump in Z only (3D fish) must NOT be filtered (filter judges
%        XY only);
%     5. pre-existing NaN frames must not be counted as "filtered";
%     6. a custom ref-point pair must drive the threshold (and the
%        reported body length), and struct arrays must be handled.

    fps = 100; nF = 400;
    mk = @() synth_fish('fps', fps, 'nFrames', nF, 'A', 0, 'U', 1, 'L', 10, ...
                        'heading', deg2rad(30));

    % ---- 1) Clean data: zero filtered ----
    fp = mk();
    out = filter_dlc_jumps(fp);
    check(out.n_points_filtered == 0, 'clean data: n_points_filtered = 0', ...
          'got %d', out.n_points_filtered);
    check_all_close(out.points, fp.points, 0, 'clean data: points untouched');
    check_close(out.jump_filter_body_length, 10, 1e-9, 'body length ref = L');

    % ---- 2) Single-frame teleport caught exactly ----
    fp2 = mk();
    fp2.points(150, 8, 1:2) = fp2.points(150, 8, 1:2) + reshape([3*fp2.gt.L, 0], 1, 1, 2);
    out2 = filter_dlc_jumps(fp2);
    check(out2.n_points_filtered == 1, 'single teleport: exactly 1 filtered', ...
          'got %d', out2.n_points_filtered);
    check(all(isnan(out2.points(150, 8, :))), 'single teleport: frame 150 point 8 NaN');
    ok_elsewhere = all(~isnan(out2.points(setdiff(1:nF,150), 8, 1:2)), 'all') && ...
                   all(~isnan(out2.points(150, setdiff(1:fp2.gt.nPts,8), 1:2)), 'all');
    check(ok_elsewhere, 'single teleport: no false positives');

    % ---- 3) 5-consecutive-frame teleport, default window 13 ----
    fp3 = mk();
    fp3.points(200:204, 12, 1) = fp3.points(200:204, 12, 1) + 3*fp3.gt.L;
    out3 = filter_dlc_jumps(fp3);
    check(out3.n_points_filtered == 5, '5-frame run: exactly 5 filtered', ...
          'got %d', out3.n_points_filtered);
    check(all(isnan(out3.points(200:204, 12, 1))), '5-frame run: all 5 bad frames NaN');
    check(~any(isnan(out3.points(199, 12, :))) && ~any(isnan(out3.points(205, 12, :))), ...
          '5-frame run: neighbors of the run untouched');
    check(all(~isnan(out3.points(setdiff(1:nF, 200:204), 12, 1))), ...
          '5-frame run: no false positives outside the run');

    % ---- 4) Z-only jump on a 3D fish is not filtered ----
    fp4 = synth_fish('fps', fps, 'nFrames', nF, 'A', 0, 'U', 1, 'L', 10, ...
                     'heading', deg2rad(30), 'z_head', 1);
    fp4.points(100, 5, 3) = fp4.points(100, 5, 3) + 50;   % huge Z teleport
    out4 = filter_dlc_jumps(fp4);
    check(out4.n_points_filtered == 0, 'Z-only jump ignored', ...
          'got %d', out4.n_points_filtered);

    % ---- 5) Pre-existing NaN frames are not counted as filtered ----
    fp5 = mk();
    fp5.points(50, 3, :) = NaN;
    out5 = filter_dlc_jumps(fp5);
    check(out5.n_points_filtered == 0, 'pre-NaN frame not counted as filtered', ...
          'got %d', out5.n_points_filtered);
    check(all(isnan(out5.points(50, 3, :))), 'pre-NaN frame stays NaN');

    % ---- 6) Custom ref pair + struct array ----
    fp6 = mk();
    fp6.points(250, 10, 1:2) = fp6.points(250, 10, 1:2) + reshape([0.3, 0], 1, 1, 2);   % 0.3 cam units
    out6 = filter_dlc_jumps([fp6, mk()], 0.5, {'p03', 'p04'});
    ref_len = fp6.gt.L / (fp6.gt.nPts - 1);      % adjacent stations: L/20 = 0.5
    check_close(out6(1).jump_filter_body_length, ref_len, 1e-9, ...
                'custom ref pair: body length = station spacing');
    check(out6(1).n_points_filtered == 1, 'custom ref pair: 0.3 jump > 0.5*ref caught', ...
          'got %d', out6(1).n_points_filtered);
    check(out6(2).n_points_filtered == 0, 'struct array: second fish clean', ...
          'got %d', out6(2).n_points_filtered);
    check(all(isnan(out6(1).points(250, 10, :))), 'custom ref pair: bad frame NaN');

    % ---- 7) Threshold is strict: a sub-threshold offset survives ----
    fp7 = mk();
    fp7.points(180, 6, 1) = fp7.points(180, 6, 1) + 0.2;   % < 0.5*L
    out7 = filter_dlc_jumps(fp7);
    check(out7.n_points_filtered == 0, 'sub-threshold offset survives', ...
          'got %d', out7.n_points_filtered);
end
