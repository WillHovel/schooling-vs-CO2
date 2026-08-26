function fish_points = filter_dlc_jumps(fish_points, jump_threshold_frac, ref_point_names, window_frames)
% FILTER_DLC_JUMPS  Remove DeepLabCut "teleport" tracking errors — a point
%                    that jumps far away for a frame or two, usually to a
%                    confusable body part or background feature — before
%                    they corrupt downstream kinematics.
%
%   fish_points = filter_dlc_jumps(fish_points, jump_threshold_frac, ref_point_names)
%   fish_points = filter_dlc_jumps(fish_points, jump_threshold_frac, ref_point_names, window_frames)
%
%   Call this right after loading (load_fish_points / load_fish_points_named),
%   BEFORE transform_fish, so the jumps don't distort the axis-fitting step.
%
%   METHOD: for each point, each frame's position is compared to the
%   MEDIAN position of a small local window of NEIGHBORING frames
%   (excluding the frame itself). If it deviates from that local median
%   by more than jump_threshold_frac * body_length, it's replaced with
%   NaN. Using a local median (not just "did it move a lot from the
%   previous frame") means a genuine fast movement across several frames
%   isn't mistaken for an error, while a single-frame teleport — which
%   stands out sharply against its own neighbors — is caught regardless
%   of whether the bad frame is the one BEFORE or AFTER a real jump.
%
%   INPUTS
%     fish_points          - struct or struct array (RAW, pre-transform_fish
%                            or post — works on whatever coordinate frame
%                            .points is in, since it only compares a point
%                            to its own recent history). Needs .points
%                            [nFrames x nPoints x nDims] and .point_names.
%     jump_threshold_frac  - (optional) fraction of body length. A point
%                            deviating from its local median by more than
%                            this fraction of body length is filtered.
%                            Default 0.5 (50%), per the "50% of the body
%                            away" heuristic.
%     ref_point_names      - (optional) {name1, name2} pair used to
%                            estimate body length (e.g. {'snout','tail'}
%                            or {'snout','peduncle'}). Default: uses the
%                            FIRST and LAST point in point_names (the
%                            usual head/tail convention elsewhere in this
%                            toolkit).
%     window_frames        - (optional) total window size (frames) used
%                            for the local median, centered on the frame
%                            being tested (excluding it). Default 13 (6
%                            frames before + 6 after) — see VALIDATION
%                            below for why this default was chosen over
%                            a smaller window. Shrinks near the start/end
%                            of the trial where a full window isn't
%                            available.
%
%   VALIDATION (tested against real DLC frogfish tracking data with
%   synthetic injected jumps, not just clean synthetic motion):
%     - Zero false positives on clean, unmodified real tracking data.
%     - A single-frame teleport is caught correctly at any window size
%       from 5 upward.
%     - A CONSECUTIVE multi-frame jump (the same point stuck on a wrong
%       feature for several frames in a row — a common real DLC failure
%       mode during brief occlusion, not just an isolated glitch) is
%       where window size matters a lot:
%         * window=5 on a 5-consecutive-frame injected jump: MISSED the
%           3 middle bad frames entirely (their neighbors were also bad,
%           so the local median agreed with them) AND false-flagged 2
%           good neighboring frames as outliers.
%         * window=9 on the same test: caught all 5 true bad frames, but
%           ALSO false-flagged the 2 good frames immediately adjacent.
%         * window=13: caught exactly the 5 true bad frames, zero false
%           positives on neighbors, zero false positives on clean data.
%     This is why the default was raised from an initial 5 to 13.
%
%   FUNDAMENTAL LIMITATION: local-median filtering can only reliably
%   catch a bad run shorter than roughly half the window's per-side
%   neighbor count. A correlated tracking failure lasting MANY consecutive
%   frames (e.g. an extended occlusion) can still evade this method,
%   because the corrupted frames start to "vote for each other" in the
%   local median. If you suspect a long bad-tracking stretch, increase
%   window_frames further, or — better — visually spot-check the point's
%   trajectory (or use DLC's own per-point likelihood column, if you kept
%   it, as an independent filter) rather than trusting this alone.
%
%   OUTPUT  fish_points — same struct, with outlier frames set to NaN in
%   .points (and .n_points_filtered / .pct_points_filtered added per fish
%   for transparency). Re-run transform_fish/compute_kinematics on the
%   result as usual — NaN handling throughout this toolkit already treats
%   missing frames correctly rather than fabricating values for them.

    if nargin < 2 || isempty(jump_threshold_frac), jump_threshold_frac = 0.5; end
    if nargin < 4 || isempty(window_frames), window_frames = 13; end
    half_w = floor(window_frames/2);

    for k = 1:numel(fish_points)
        pts = fish_points(k).points;   % [nFrames x nPoints x nDims]
        [nFrames, nPoints, nDims] = size(pts);
        pn = fish_points(k).point_names;

        % ---- Body length reference ----
        if nargin >= 3 && ~isempty(ref_point_names)
            i1 = find(strcmpi(pn, ref_point_names{1}), 1);
            i2 = find(strcmpi(pn, ref_point_names{2}), 1);
        else
            i1 = 1; i2 = nPoints;
        end
        if isempty(i1) || isempty(i2)
            warning(['filter_dlc_jumps: reference points not found for fish %d ("%s") — ' ...
                     'skipping (no filtering applied to this fish).'], k, safe_name(fish_points(k)));
            continue;
        end

        p1 = squeeze(pts(:, i1, 1:min(2,nDims)));
        p2 = squeeze(pts(:, i2, 1:min(2,nDims)));
        body_len_per_frame = sqrt(sum((p1-p2).^2, 2));
        body_length = median(body_len_per_frame, 'omitnan');

        if isnan(body_length) || body_length <= 0
            warning(['filter_dlc_jumps: could not estimate a valid body length for fish %d ' ...
                     '("%s") — skipping (no filtering applied to this fish).'], ...
                     k, safe_name(fish_points(k)));
            continue;
        end
        threshold = jump_threshold_frac * body_length;

        % ---- Per-point local-median outlier filtering ----
        % Tracked point-frames (XY present) — excludes leading/trailing
        % all-empty padding rows common in exported CSVs.
        n_present = sum(any(~isnan(pts(:, :, 1:min(2,nDims))), 3), 'all');
        n_filtered = 0;
        for pi = 1:nPoints
            for f = 1:nFrames
                lo = max(1, f-half_w);
                hi = min(nFrames, f+half_w);
                window_idx = [lo:f-1, f+1:hi];
                if numel(window_idx) < 2, continue; end   % not enough neighbors to judge

                window_pos = squeeze(pts(window_idx, pi, 1:min(2,nDims)));
                if size(window_pos,1) == 1, window_pos = window_pos'; end
                local_median = median(window_pos, 1, 'omitnan');
                if any(isnan(local_median)), continue; end

                cur_pos = squeeze(pts(f, pi, 1:min(2,nDims)))';
                if any(isnan(cur_pos)), continue; end

                dev = norm(cur_pos - local_median);
                if dev > threshold
                    pts(f, pi, :) = NaN;
                    n_filtered = n_filtered + 1;
                end
            end
        end

        fish_points(k).points = pts;
        fish_points(k).n_points_filtered  = n_filtered;
        fish_points(k).pct_points_filtered = 100 * n_filtered / max(n_present, 1);
        fish_points(k).jump_filter_body_length  = body_length;
        fish_points(k).jump_filter_threshold_frac = jump_threshold_frac;

        fprintf(['filter_dlc_jumps: %s | body length ref = %.4g | threshold = %.4g (%.0f%% of body) | ' ...
                 '%d/%d tracked point-frames filtered (%.2f%%)\n'], ...
                safe_name(fish_points(k)), body_length, threshold, 100*jump_threshold_frac, ...
                n_filtered, n_present, fish_points(k).pct_points_filtered);
    end
end

function s = safe_name(fp)
    if isfield(fp, 'name') && ~isempty(fp.name), s = fp.name; else, s = 'unnamed'; end
end
