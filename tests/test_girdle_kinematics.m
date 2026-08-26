function test_girdle_kinematics()
% TEST_GIRDLE_KINEMATICS  Ground-truth checks for compute_girdle_kinematics.m.
%
%   A girdle point rigidly attached to a straight (A=0), heading-0 fish at
%   body-frame position (x_b(t) = -3 + 0.1*sin(2*pi*f0*t), y_b = 0.5) — a
%   fore-aft oscillation of 0.1 camera units on a 10-unit body. After
%   projection into the body frame (head at X=0, tail at X=1, X = -x_b/L):
%     X = 0.3 - 0.01*sin(2*pi*f0*t),  Y = 0.05 exactly
%   so protraction_range = 0.02 BL, lateral_range = 0, and the girdle
%   frequency is exactly f0 = 2 Hz (FFT bin exact at nF=400, fps=100).
%   Also checks the error paths: frame-count mismatch, missing column,
%   multi-animal struct array.

    d = tempname;
    mkdir(d);
    cleanup = onCleanup(@() rmdir(d, 's'));
    f = @(nm) fullfile(d, nm);

    fps = 100; nF = 400; f0 = 2;
    fp = synth_fish('fps', fps, 'nFrames', nF, 'A', 0, 'U', 1, 'L', 10, 'heading', 0);
    tf = transform_fish(fp);
    t = fp.gt.t;

    % girdle point rigidly attached in BODY coords, re-expressed in camera
    % coords exactly as synth_fish puts the fish on camera:
    %   cam = center + U*[cos h0; sin h0]*t + R(h)*(x_b, y_b)
    xb = -3 + 0.1*sin(2*pi*f0*t);
    yb = 0.5;
    hd = fp.gt.heading_fn(t);
    gx = fp.gt.center(1) + fp.gt.U*cos(fp.gt.heading_fn(0))*t + cos(hd).*xb - sin(hd).*yb;
    gy = fp.gt.center(2) + fp.gt.U*sin(fp.gt.heading_fn(0))*t + sin(hd).*xb + cos(hd).*yb;
    gz = zeros(nF, 1);
    rows = [ {'GP_X','GP_Y','GP_Z'}; num2cell([gx, gy, gz]) ];
    write_csv(f('girdle.csv'), rows);

    g = compute_girdle_kinematics(f('girdle.csv'), 'GP', tf, fps, 0.5);
    check_all_close(g.X, 0.3 - 0.01*sin(2*pi*f0*t), 1e-9, 'X = body-frame fore-aft');
    check_all_close(g.Y, 0.05*ones(nF,1), 1e-9, 'Y = 0.5/L = 0.05 BL');
    % the sampled grid misses the exact oscillation peak: the measured
    % range is 0.02 * max(sin(2*pi*f0*t)) on t = (0:nF-1)/fps
    check_close(g.protraction_range_BL, 0.02 * max(sin(2*pi*f0*t)), 1e-9, ...
                'protraction range = 0.02 BL');
    check_close(g.lateral_range_BL, 0, 1e-9, 'lateral range = 0 (fixed y)');
    check_close(g.girdle_freq_Hz, f0, 1e-9, 'girdle freq = 2 Hz');
    check(g.n_valid == nF, 'all frames valid', 'got %d', g.n_valid);
    check(g.pct_valid == 100, 'pct_valid = 100', 'got %.1f', g.pct_valid);

    % ---- Frame-count mismatch ----
    tf2 = tf;
    tf2.transform_params = tf.transform_params(1:end-1);
    threw = false;
    try
        compute_girdle_kinematics(f('girdle.csv'), 'GP', tf2, fps, 0.5);
    catch err
        threw = true;
        check(contains(err.message, 'frame count mismatch'), 'mismatch error', ...
              'msg = %s', err.message);
    end
    check(threw, 'frame-count mismatch must error');

    % ---- Missing column ----
    threw = false;
    try
        compute_girdle_kinematics(f('girdle.csv'), 'Nope', tf, fps, 0.5);
    catch err
        threw = true;
        check(contains(err.message, 'not found'), 'missing column error', ...
              'msg = %s', err.message);
    end
    check(threw, 'missing column must error');

    % ---- Multi-animal struct array ----
    threw = false;
    try
        compute_girdle_kinematics(f('girdle.csv'), 'GP', [tf tf], fps, 0.5);
    catch err
        threw = true;
        check(contains(err.message, 'single-animal'), 'multi-animal error', ...
              'msg = %s', err.message);
    end
    check(threw, 'struct array must error');
end
