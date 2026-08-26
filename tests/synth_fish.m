function fp = synth_fish(varargin)
% SYNTH_FISH  Build a synthetic fish in CAMERA coordinates with exactly
%             known kinematics, for ground-truth testing.
%
%   fp = synth_fish('Name', Value, ...)
%
%   The fish is a straight spine of body length L with nPts stations
%   s_i = linspace(0,1,nPts) (head at s=0). In the BODY frame the spine
%   lies on the negative x-axis (head at the origin, tail at -L) and each
%   station carries a lateral traveling wave
%       y_b(s,t) = A(s) * sin(2*pi*(f0*t - s/lambda))
%   (phase lags toward the tail — the propulsive direction). A(s) may be
%   a scalar (constant amplitude — the exactly-validated case) or a
%   function handle of s. The fish swims along a path at ground speed U:
%     'path' = 'straight'  -> heading(t) fixed, position = U*t along heading
%              'stationary'-> heading(t) may still vary, no translation
%              'circle'    -> heading(t) = heading0 + (U/R)*t, position on a
%                             circle of radius R centered at 'center'
%   Camera coords: p_i(t) = c(t) + R(heading(t)) * [x_b, y_b]  (2D), with an
%   optional vertical (z) term:  z_i(t) = z_head(t) + z_off_i, where
%   z_head is a scalar height or function of t and z_off_i is a per-point
%   offset (e.g. a pitch ramp). Any z option makes the fish 3D (has_z).
%
%   OUTPUT fp — struct in the toolkit's fish_points schema:
%     .name, .frames, .point_names, .points, .has_z, .format
%   plus .gt (ground truth): f0, lambda, A (fn of s), U, L, heading (fn of
%   t), path, nPts — everything a test needs to compute expected metrics.
%
%   Name/Value options:
%     fps (100)  nFrames (400)  nPts (21)
%     f0 (2 Hz)          beat frequency of the traveling wave
%     lambda (1 BL)      propulsive wavelength in BL
%     A (0.05)           wave amplitude, in CAMERA units, scalar or
%                        function handle of s
%     U (1)              ground speed in camera units per second
%     L (10)             body length in camera units
%     heading (0)        scalar (rad) or function handle of t (rad)
%     path ('straight')  'straight' | 'stationary' | 'circle'
%     center ([0 0])     circle center (camera units)
%     z_head ([])        scalar height or function handle of t -> enables 3D
%     z_off ([])         [1 x nPts] per-point vertical offsets (camera units)
%     jitter (0)         gaussian tracking noise, camera units (std)
%
%   NOTE: amplitude A is in camera units; the toolkit's metrics are in BL,
%   so expected amplitudes are A/L (tests divide accordingly).

    p = inputParser;
    addParameter(p, 'fps',      100);
    addParameter(p, 'nFrames',  400);
    addParameter(p, 'nPts',     21);
    addParameter(p, 'f0',       2);
    addParameter(p, 'lambda',   1);
    addParameter(p, 'A',        0.05);
    addParameter(p, 'U',        1);
    addParameter(p, 'L',        10);
    addParameter(p, 'heading',  0);
    addParameter(p, 'path',     'straight');
    addParameter(p, 'center',   [0 0]);
    addParameter(p, 'z_head',   []);
    addParameter(p, 'z_off',    []);
    addParameter(p, 'jitter',   0);
    parse(p, varargin{:});
    o = p.Results;

    t    = (0:o.nFrames-1)' / o.fps;
    s    = linspace(0, 1, o.nPts);

    A_fn = o.A;
    if isnumeric(A_fn), A_fn = @(ss) o.A + 0*ss; end

    % heading(t) and path position c(t) in camera coords
    if isnumeric(o.heading)
        heading_fn = @(tt) o.heading + 0*tt;
    else
        heading_fn = o.heading;
    end
    switch lower(o.path)
        case 'straight'
            h0 = heading_fn(0);
            c_fn = @(tt) o.center(:) + o.U .* [cos(h0); sin(h0)] .* tt(:)';
            heading_fn = @(tt) h0 + 0*tt(:);           % constant heading
        case 'stationary'
            c_fn = @(tt) repmat(o.center(:), 1, numel(tt));
        case 'circle'
            % c(t) = center + R*[sin h, -cos h]  ->  c' = R*h'*[cos h, sin h],
            % so R = U/omega makes the ground speed exactly U along heading.
            h0 = heading_fn(0); h1 = heading_fn(1e-3);
            omega = (h1 - h0) / 1e-3;
            if abs(omega) < eps
                error('synth_fish: circle path needs a time-varying heading function handle.');
            end
            R_c = o.U / omega;
            % cstack keeps the [sin; -cos] stack 2 x numel(tt) for vector
            % tt (a plain [col; col] would stack into one 2*numel column)
            c_fn = @(tt) o.center(:) + R_c .* cstack(heading_fn(tt));
        otherwise
            error('synth_fish: unknown path "%s"', o.path);
    end

    has_z = ~isempty(o.z_head) || ~isempty(o.z_off);
    nDims = 2 + has_z;
    pts = NaN(o.nFrames, o.nPts, nDims);

    for i = 1:o.nPts
        x_b = -s(i) * o.L;                              % head at 0, tail at -L
        y_b = A_fn(s(i)) * sin(2*pi*(o.f0 .* t - s(i)/o.lambda));
        hd  = heading_fn(t);
        cx  = c_fn(t);
        pts(:, i, 1) = cx(1,:)' + cos(hd) .* x_b - sin(hd) .* y_b;
        pts(:, i, 2) = cx(2,:)' + sin(hd) .* x_b + cos(hd) .* y_b;
        if has_z
            if isnumeric(o.z_head)
                zt = o.z_head + 0*t;
            else
                zt = o.z_head(t);
            end
            if isempty(o.z_off)
                pts(:, i, 3) = zt;
            else
                pts(:, i, 3) = zt + o.z_off(i);
            end
        end
    end

    if o.jitter > 0
        rng(1, 'twister');
        pts = pts + o.jitter * randn(size(pts));
    end

    fp.name        = 'synthetic';
    fp.frames      = (1:o.nFrames)';
    fp.point_names = arrayfun(@(k) sprintf('p%02d', k), 1:o.nPts, 'UniformOutput', false);
    fp.points      = pts;
    fp.has_z       = has_z;
    fp.format      = 'synthetic';

    fp.gt          = o;
    fp.gt.A_fn     = A_fn;
    fp.gt.heading_fn = heading_fn;
    fp.gt.s        = s;
    fp.gt.t        = t;
end

function uv = cstack(hd)
% Unit circle direction vectors as a 2 x numel(hd) matrix.
    hd = hd(:).';
    uv = [sin(hd); -cos(hd)];
end
