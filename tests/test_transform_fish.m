function test_transform_fish()
% TEST_TRANSFORM_FISH  Ground-truth checks for transform_fish.m.
%
% Ground truth strategy: synth_fish builds fish with EXACTLY known camera
% geometry, so every property of the transform is predictable:
%   1. camera roundtrip — X/Y + transform_params must reconstruct the raw
%      camera coordinates to floating-point precision (this pins down the
%      whole rotation/translation/scaling/mirroring chain at once);
%   2. a perfectly straight fish must land on Y = 0 with X = station
%      position (0=head, 1=tail);
%   3. rigid rotations/translations/scaling of the camera must leave the
%      normalized X/Y output unchanged (frame-invariance);
%   4. reversed point order must be auto-mirrored (point_names{1} stays
%      the head, X=0..1, lateral excursions preserved);
%   5. bl_override must replace the measured body length;
%   6. missing endpoints must leave the frame NaN and be counted in
%      pct_frames_valid;
%   7. pre_transformed (CURVES) data must pass through untouched;
%   8. a 3-point fish (single middle point) must fall back to the
%      head-to-tail chord fit.

    fps = 100; nF = 400;

    % ---- 1) Camera roundtrip on a swimming, rotated, wavy fish ----
    fp = synth_fish('fps', fps, 'nFrames', nF, 'f0', 2, 'lambda', 1, ...
                    'A', 0.05, 'U', 1, 'L', 10, 'heading', deg2rad(30));
    tf = transform_fish(fp);
    check(tf.pct_frames_valid == 100, 'pct_frames_valid', 'got %.2f', tf.pct_frames_valid);
    % synth_fish lays the body along the NEGATIVE x-axis (tail at smaller
    % rotated x than the head), so the 0->1 mirroring fires every frame
    check(tf.n_frames_reversed == nF, 'head-first synth fish mirrors every frame');

    tp = tf.transform_params;
    maxerr = 0;
    for f = 1:nF
        th = tp(f).theta; a = tp(f).a; bl = tp(f).bl; x1 = tp(f).x1; sf = tp(f).sign_flip;
        for i = 1:size(tf.X, 2)
            xcam = (sf .* tf.X(f,i) .* bl + x1) .* cos(th) + tf.Y(f,i) .* bl .* sin(th);
            ycam = (tf.Y(f,i) .* bl - xcam .* sin(th)) / cos(th) + a;
            maxerr = max(maxerr, abs(xcam - fp.points(f,i,1)));
            maxerr = max(maxerr, abs(ycam - fp.points(f,i,2)));
        end
    end
    check(maxerr < 1e-9, 'camera roundtrip error', 'got %.3e', maxerr);
    check_all_close(tf.X(:,1),    zeros(nF,1), 1e-9, 'X(head) = 0');
    check_all_close(tf.X(:,end),  ones(nF,1),  1e-9, 'X(tail) = 1');
    % bl_per_frame comes from the LS line fit over the wavy midline, so it
    % tracks L only to fit precision (~1e-4); the straight fish below is exact
    check_all_close(tf.bl_per_frame, fp.gt.L*ones(nF,1), 1e-3, 'bl_per_frame ~ true body length');

    % ---- 2) Straight fish: Y = 0, X = station position ----
    fs2 = synth_fish('fps', fps, 'nFrames', nF, 'A', 0, 'U', 1, 'L', 10, ...
                     'heading', deg2rad(-40), 'center', [3; 7]);
    tf2 = transform_fish(fs2);
    check_all_close(tf2.Y, zeros(nF, size(tf2.Y,2)), 1e-9, 'straight fish Y = 0');
    check_all_close(tf2.X, repmat(linspace(0,1,fs2.gt.nPts), nF, 1), 1e-9, 'straight fish X = s');
    check_all_close(tf2.bl_per_frame, fs2.gt.L*ones(nF,1), 1e-9, 'straight fish bl_per_frame = L exactly');

    % ---- 3) Rigid camera motion + scaling invariance ----
    fps3 = synth_fish('fps', fps, 'nFrames', nF, 'f0', 2, 'lambda', 1, 'A', 0.05, ...
                      'U', 1, 'L', 10, 'heading', deg2rad(30));
    cam_rot = [cosd(50) -sind(50); sind(50) cosd(50)];
    p2 = reshape(fps3.points(:, :, 1:2), [], 2);
    p2 = 1.7 * (p2 * cam_rot' + [3 7]);
    fps3.points(:, :, 1:2) = reshape(p2, nF, [], 2);
    tf3 = transform_fish(fps3);
    % exact frame-invariance would hold in real arithmetic; the LS midline
    % fit is only rigid-invariant to ~1e-5 in X on this wavy fish. Y is
    % worse (~4e-4): the fit minimizes VERTICAL residuals in camera coords,
    % which is not rotation-invariant on a wavy midline.
    check_all_close(tf3.X, tf.X, 1e-5, 'X invariant under rigid camera motion + scale');
    check_all_close(tf3.Y, tf.Y, 1e-3, 'Y invariant under rigid camera motion + scale');

    % ---- 4) Reversed point order -> auto-mirror, head/tail preserved ----
    fps4 = synth_fish('fps', fps, 'nFrames', nF, 'f0', 2, 'lambda', 1, 'A', 0.05, ...
                      'U', 1, 'L', 10, 'heading', deg2rad(10));
    fps4.points    = flip(fps4.points, 2);
    fps4.point_names = fliplr(fps4.point_names);
    tf4 = transform_fish(fps4);
    check(tf4.n_frames_reversed == 0, 'reversed order needs no mirroring', ...
          'got %d/%d', tf4.n_frames_reversed, nF);
    check_all_close(tf4.X(:,1),   zeros(nF,1), 1e-9, 'reversed fish: X(head) = 0');
    check_all_close(tf4.X(:,end), ones(nF,1),  1e-9, 'reversed fish: X(tail) = 1');
    % Same physical point (old tail = new head) must have the same lateral
    % trace. Only to ~1e-4: the two fish sit at different camera headings
    % (10 vs 30 deg), and the LS midline fit's vertical residuals are not
    % rigid-invariant on a wavy midline (same effect as check 3 above).
    check_all_close(tf4.Y(:,1), tf.Y(:,end), 1e-4, 'mirroring preserves lateral excursions');

    % ---- 5) bl_override ----
    tf5 = transform_fish(fp, 5);
    check_all_close(tf5.bl_per_frame, 5*ones(nF,1), 1e-9, 'bl_override applied');
    % X(tail) = measured chord length / override; the chord length on the
    % wavy fish is only good to ~1e-4 (LS midline fit, same as check 1)
    check_all_close(tf5.X(:,end), 2*ones(nF,1), 1e-3, 'X(tail) = true_length/override');

    % ---- 6) Missing endpoints -> NaN frames, correct validity % ----
    fps6 = fp;
    fps6.points(10:12, 1, :) = NaN;      % head missing on 3 frames
    fps6.points(20:22, end, :) = NaN;    % tail missing on 3 frames
    tf6 = transform_fish(fps6);
    check_all_close(tf6.X(10:12, :), NaN(3, size(tf6.X,2)), 1, 'head-missing frames stay NaN');
    check_all_close(tf6.X(20:22, :), NaN(3, size(tf6.X,2)), 1, 'tail-missing frames stay NaN');
    check_close(tf6.pct_frames_valid, 100*(nF-6)/nF, 1e-9, 'pct_frames_valid counts NaN frames');

    % ---- 7) pre_transformed passthrough (CURVES) ----
    fpp = fp;
    fpp.pre_transformed = true;
    fpp.X = rand(nF, fp.gt.nPts);
    fpp.Y = rand(nF, fp.gt.nPts);
    fpp_saved = fpp;
    tpp = transform_fish(fpp);
    check(isequaln(tpp, fpp_saved), 'pre_transformed data passes through unchanged');

    % ---- 8) 3-point fish falls back to chord fit ----
    fp8 = synth_fish('fps', fps, 'nFrames', nF, 'nPts', 3, 'A', 0, 'U', 1, 'L', 10, ...
                     'heading', deg2rad(25));
    tf8 = transform_fish(fp8);
    check_all_close(tf8.Y, zeros(nF,3), 1e-9, '3-point straight fish Y = 0');
    check_all_close(tf8.X, repmat([0 0.5 1], nF, 1), 1e-9, '3-point straight fish X = 0/0.5/1');
end
