function test_apply_body_transform()
% TEST_APPLY_BODY_TRANSFORM  Ground-truth checks for apply_body_transform.m
%   (the projection of arbitrary points into transform_fish's body frame).
%
%   1. Self-consistency: projecting each MIDLINE point must reproduce the
%      exact X/Y that transform_fish computed for it (same formula, so
%      this pins the interface, any drift between the two breaks girdle
%      and fin analyses).
%   2. A point rigidly attached to a straight fish has CONSTANT body-frame
%      coordinates equal to its known body position (head at X=0).
%   3. Frames where the midline transform failed must be NaN.

    fps = 100; nF = 400;
    fp = synth_fish('fps', fps, 'nFrames', nF, 'f0', 2, 'lambda', 1, 'A', 0.05, ...
                    'U', 1, 'L', 10, 'heading', deg2rad(30));
    tf = transform_fish(fp);
    tp = tf.transform_params;

    % ---- 1) Self-consistency on every midline point (one point at a time,
    %         apply_body_transform's interface is [nFrames x 2/3]) ----
    for i = 1:size(tf.X, 2)
        [Xb, Yb] = apply_body_transform(tf.points(:, i, 1:2), tp);
        check_all_close(Xb, tf.X(:, i), 1e-9, 'apply_body_transform X == transform_fish X');
        check_all_close(Yb, tf.Y(:, i), 1e-9, 'apply_body_transform Y == transform_fish Y');
    end

    % ---- 2) Rigid fin-base point: constant, known body coordinates ----
    fs = synth_fish('fps', fps, 'nFrames', nF, 'A', 0, 'U', 1, 'L', 10, ...
                    'heading', deg2rad(20));
    tfs = transform_fish(fs);
    bx = -3; by = 0.5;                       % body coords of the fin base
    hd = fs.gt.heading_fn(fs.gt.t);
    cx = fs.gt.center(1) + fs.gt.U * cos(fs.gt.heading_fn(0)) * fs.gt.t;
    cy = fs.gt.center(2) + fs.gt.U * sin(fs.gt.heading_fn(0)) * fs.gt.t;
    fin_raw = [cx + cos(hd)*bx - sin(hd)*by, cy + sin(hd)*bx + cos(hd)*by];
    [Xf, Yf] = apply_body_transform(fin_raw, tfs.transform_params);
    check_all_close(Xf, (-bx/fs.gt.L)*ones(nF,1), 1e-9, 'rigid fin point X = -bx/L');
    check_all_close(Yf, ( by/fs.gt.L)*ones(nF,1), 1e-9, 'rigid fin point Y = by/L');

    % ---- 3) Failed midline frames -> NaN ----
    fp3 = fp;
    fp3.points(50:52, 1, :) = NaN;           % kills the midline fit there
    tf3 = transform_fish(fp3);
    [X3, Y3] = apply_body_transform(tf3.points(:, :, 1:2), tf3.transform_params);
    check(all(isnan(X3(50:52))) && all(isnan(Y3(50:52))), ...
          'failed midline frames project to NaN');
    check(~any(isnan(X3(60:70))), 'valid frames still project');
end
