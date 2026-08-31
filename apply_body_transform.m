function [X, Y, Z] = apply_body_transform(raw_xyz, transform_params)
% APPLY_BODY_TRANSFORM  Project an arbitrary tracked point into the SAME
%                       body-relative (0=head->1=tail, BL-normalized) frame
%                       that transform_fish computed for the midline.
%
%   [X, Y, Z] = apply_body_transform(raw_xyz, transform_params)
%
%   Use this to analyze a point that ISN'T one of the midline points fed
%   into transform_fish, e.g. a fin root/girdle marker, in body-relative
%   coordinates, so its motion is comparable frame-to-frame regardless of
%   the fish's own translation/rotation/size. This is what makes girdle
%   protraction-retraction analysis possible (compute_girdle_kinematics.m).
%
%   INPUTS
%     raw_xyz          - [nFrames x 2] or [nFrames x 3] raw (untransformed)
%                         coordinates of the point of interest, same units
%                         and frame indexing as what was passed to
%                         transform_fish.
%     transform_params - fish_points(i).transform_params from transform_fish
%                         (a 1 x nFrames struct array with theta, a, x1,
%                         bl, sign_flip per frame, NaN where that frame's
%                         midline transform failed).
%
%   OUTPUTS
%     X, Y, Z  [nFrames x 1]  body-relative position of the point, in the
%              same 0(head)->1(tail) BL-normalized convention as the
%              midline's own .X/.Y/.Z. Z is only returned if raw_xyz has
%              3 columns; otherwise it's returned empty.
%
%   NOTE: a frame is NaN here whenever the MIDLINE's own transform failed
%   for that frame (transform_params(f).bl is NaN) OR this point itself is
%   NaN in that frame: you can't express a point in a coordinate frame
%   that itself couldn't be computed.

    nFrames = size(raw_xyz, 1);
    has_z   = size(raw_xyz, 2) >= 3;

    X = NaN(nFrames, 1);
    Y = NaN(nFrames, 1);
    Z = [];
    if has_z, Z = NaN(nFrames, 1); end

    for f = 1:nFrames
        tp = transform_params(f);
        % isempty() catches the case where this frame's transform_params
        % element was never explicitly assigned (see transform_fish.m's
        % CHANGE NOTE: struct-array preallocation quirks can leave
        % unset elements as [] rather than NaN). Checking both makes this
        % robust even if that happens again for some other reason.
        if isempty(tp.bl) || isnan(tp.bl), continue; end

        x = raw_xyz(f,1); y = raw_xyz(f,2);
        if isnan(x) || isnan(y), continue; end

        x_r = x*cos(tp.theta) - (y - tp.a)*sin(tp.theta);
        y_r = x*sin(tp.theta) + (y - tp.a)*cos(tp.theta) + tp.a;

        X(f) = tp.sign_flip * (x_r - tp.x1) / tp.bl;
        Y(f) = (y_r - tp.a) / tp.bl;

        if has_z
            z = raw_xyz(f,3);
            if ~isnan(z)
                % Z is expressed relative to this point's OWN first valid
                % frame (there's no single universal DV reference point
                % once you're off the midline), see compute_girdle_kinematics
                % for how this is used in practice.
                Z(f) = z / tp.bl;
            end
        end
    end
end
