function stance = compute_stance_swing(fin, varargin)
% COMPUTE_STANCE_SWING  Estimate contact time / duty factor from fin-tip
%                       velocity — for WALKING/PUNTING gaits (e.g.
%                       Polypterus, epaulette shark "punting"), not
%                       free-swimming sculling.
%
%   stance = compute_stance_swing(fin)
%   stance = compute_stance_swing(fin, 'Name', Value, ...)
%
%   METHOD (IMPORTANT — READ BEFORE TRUSTING THIS)
%   There's no force/pressure sensor here, so "contact" is a KINEMATIC
%   PROXY: a frame is classified as STANCE (fin planted) when the tip's
%   speed drops below a threshold fraction of its own 95th-percentile
%   speed for this trial, and SWING (repositioning) otherwise. This is a
%   standard proxy in locomotion kinematics, but it is NOT the same
%   measurement as a true contact sensor / pressure-mat duty factor. If
%   you have real contact-time data (as Shoval mentioned), TREAT THAT AS
%   GROUND TRUTH and use this as a secondary/exploratory check, not a
%   replacement — the two can disagree, especially if the fin doesn't
%   fully stop during stance (common in a slow punt vs a hard plant).
%
%   INPUTS
%     fin       - output struct from compute_fin_kinematics.m (needs
%                 .tip_speed, .fps, .valid).
%     Name/Value options:
%       'ThresholdFrac'   - stance = tip_speed < ThresholdFrac * p95(tip_speed).
%                            Default 0.25. Lower = stricter (less generous
%                            about calling something "stance").
%       'MinBoutFrames'   - discard bouts shorter than this many frames
%                            (noise/jitter filter). Default 3.
%       'MaxGapFrames'    - a data dropout (NaN) shorter than this many
%                            frames is bridged (treated as continuing
%                            whatever phase surrounded it) rather than
%                            splitting a bout. Longer gaps end the bout
%                            and are excluded. Default 2.
%
%   OUTPUT  stance — struct:
%     .is_stance          [nFrames x 1] logical (NaN frames -> false)
%     .threshold_used      the actual speed threshold (same units as tip_speed)
%     .stance_bouts_frames  [nBouts x 2] start/end frame indices
%     .swing_bouts_frames   [nBouts x 2] start/end frame indices
%     .mean_contact_time_s  mean stance bout duration (s)
%     .mean_swing_time_s    mean swing bout duration (s)
%     .duty_factor          mean_contact_time / (mean_contact_time + mean_swing_time)
%     .n_cycles             number of complete stance->swing->stance cycles found
%     .pct_valid_consecutive  % of frames with usable frame-to-frame speed —
%                              LOW VALUES HERE (<70%) MEAN THIS ESTIMATE IS
%                              LIKELY UNRELIABLE due to tracking gaps
%                              fragmenting real bouts. Check this before
%                              trusting duty_factor.

    p = inputParser;
    addParameter(p, 'ThresholdFrac', 0.25);
    addParameter(p, 'MinBoutFrames', 3);
    addParameter(p, 'MaxGapFrames', 2);
    parse(p, varargin{:});
    thresh_frac    = p.Results.ThresholdFrac;
    min_bout       = p.Results.MinBoutFrames;
    max_gap        = p.Results.MaxGapFrames;

    tip_speed = fin.tip_speed;
    fps       = fin.fps;
    nFrames   = numel(tip_speed);

    % tip_speed is 0 (not NaN) for frames where step_dist couldn't be
    % computed (see compute_fin_kinematics.m — step_dist initialized to
    % zeros, only filled where BOTH this and prior frame are valid).
    % Treat those zero-by-construction frames as "unknown", not "stance",
    % or every tracking gap would masquerade as a plant.
    known = fin.valid & [false; fin.valid(1:end-1)];   % both this & prior frame tracked
    speed_known = tip_speed;
    speed_known(~known) = NaN;

    pct_valid_consecutive = 100 * sum(known) / nFrames;
    if pct_valid_consecutive < 70
        warning(['compute_stance_swing: only %.1f%% of frames have usable ' ...
                 'frame-to-frame tip speed (rest are tracking gaps). Stance/swing ' ...
                 'bouts are likely fragmented by missing data, not real transitions — ' ...
                 'treat duty_factor as a rough estimate, not a precise measurement, ' ...
                 'until tracking coverage improves.'], pct_valid_consecutive);
    end

    if sum(known) < 10
        warning('compute_stance_swing: too few valid frames (%d) to estimate duty factor reliably.', sum(known));
        stance = empty_result(nFrames, NaN, pct_valid_consecutive);
        return;
    end

    p95 = prctile(speed_known(known), 95);
    threshold = thresh_frac * p95;

    is_stance = false(nFrames, 1);
    is_stance(known) = speed_known(known) < threshold;

    % ---- Bridge short data gaps (<= max_gap frames) rather than splitting bouts ----
    phase = NaN(nFrames,1);   % 1 = stance, 0 = swing, NaN = unknown/gap
    phase(known) = is_stance(known);

    i = 1;
    while i <= nFrames
        if isnan(phase(i))
            j = i;
            while j <= nFrames && isnan(phase(j)), j = j + 1; end
            gap_len = j - i;
            if gap_len <= max_gap && i > 1 && j <= nFrames && phase(i-1) == phase(j)
                phase(i:j-1) = phase(i-1);   % bridge — same phase on both sides
            end
            i = j;
        else
            i = i + 1;
        end
    end

    is_stance_bridged = (phase == 1);
    is_swing_bridged  = (phase == 0);

    [stance_bouts, stance_durs] = find_bouts(is_stance_bridged, min_bout, fps);
    [swing_bouts,  swing_durs]  = find_bouts(is_swing_bridged,  min_bout, fps);

    if isempty(stance_durs) || isempty(swing_durs)
        warning(['compute_stance_swing: no complete stance and/or swing bouts survived ' ...
                 'the MinBoutFrames filter (%d frames). Try lowering MinBoutFrames or ' ...
                 'check that this trial actually has a walking/punting gait (this method ' ...
                 'assumes distinct plant/lift phases — it is not meaningful for continuous ' ...
                 'sculling or free-swimming fin motion).'], min_bout);
        mean_contact = NaN; mean_swing = NaN; duty_factor = NaN;
    else
        mean_contact = mean(stance_durs);
        mean_swing   = mean(swing_durs);
        duty_factor  = mean_contact / (mean_contact + mean_swing);
    end

    n_cycles = min(size(stance_bouts,1), size(swing_bouts,1));

    stance.is_stance            = is_stance_bridged;
    stance.threshold_used        = threshold;
    stance.stance_bouts_frames   = stance_bouts;
    stance.swing_bouts_frames    = swing_bouts;
    stance.mean_contact_time_s   = mean_contact;
    stance.mean_swing_time_s     = mean_swing;
    stance.duty_factor           = duty_factor;
    stance.n_cycles              = n_cycles;
    stance.pct_valid_consecutive = pct_valid_consecutive;

    fprintf(['compute_stance_swing: %s->%s | %.1f%% frames usable | %d cycles | ' ...
             'contact=%s s  swing=%s s  duty factor=%s\n'], ...
            fin.rootName, fin.tipName, pct_valid_consecutive, n_cycles, ...
            fmt_val(mean_contact), fmt_val(mean_swing), fmt_val(duty_factor));
    if pct_valid_consecutive < 70
        fprintf('  ^ LOW COVERAGE WARNING ABOVE APPLIES — verify against your measured contact-time CSV.\n');
    end
end


function [bouts, durs_s] = find_bouts(mask, min_bout_frames, fps)
    bouts = [];
    durs_s = [];
    i = 1; n = numel(mask);
    while i <= n
        if mask(i)
            j = i;
            while j <= n && mask(j), j = j + 1; end
            bout_len = j - i;
            if bout_len >= min_bout_frames
                bouts = [bouts; i, j-1]; %#ok<AGROW>
                durs_s = [durs_s; bout_len / fps]; %#ok<AGROW>
            end
            i = j;
        else
            i = i + 1;
        end
    end
end

function stance = empty_result(nFrames, thresh, pct_valid)
    stance.is_stance = false(nFrames,1);
    stance.threshold_used = thresh;
    stance.stance_bouts_frames = [];
    stance.swing_bouts_frames = [];
    stance.mean_contact_time_s = NaN;
    stance.mean_swing_time_s = NaN;
    stance.duty_factor = NaN;
    stance.n_cycles = 0;
    stance.pct_valid_consecutive = pct_valid;
end

function s = fmt_val(v)
    if isnan(v), s = 'NaN'; else, s = sprintf('%.4f', v); end
end
