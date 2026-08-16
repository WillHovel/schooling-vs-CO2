function cyc = compute_cycle_frequencies(y, fps, min_freq, max_freq, ref_freq_Hz)
% COMPUTE_CYCLE_FREQUENCIES  Per-cycle tail-beat frequency, not a single
%                             trial-averaged Hz value.
%
%   cyc = compute_cycle_frequencies(y, fps)
%   cyc = compute_cycle_frequencies(y, fps, min_freq, max_freq)
%   cyc = compute_cycle_frequencies(y, fps, min_freq, max_freq, ref_freq_Hz)
%
%   compute_kinematics.m's head_TBF/tail_TBF give ONE FFT-averaged
%   frequency for the whole trial. This instead detects each individual
%   oscillation cycle via interpolated positive-going zero-crossings and
%   reports one frequency PER CYCLE — e.g. 4 beats in a clip -> 4 numbers,
%   not 1 average.
%
%   IMPORTANT ACCURACY NOTE: body-relative Y time series in this toolkit
%   often carry slow drift/trend (e.g. from imperfect frame-by-frame axis
%   fitting) rather than being a clean signal oscillating around zero.
%   Naive zero-crossing counting on such a signal badly under/over-counts
%   cycles, so this function high-pass detrends first (subtracts a moving-
%   average baseline). Even so, on at least one validated trial in this
%   toolkit's test data, cycle-based detection returned a NOTICEABLY
%   different frequency (~0.21 Hz, 1 cycle found) than the independently
%   validated FFT-based value for the same signal (~0.30 Hz, confirmed
%   against both a full power-spectrum check AND a visual beat count).
%   This is a genuine accuracy limitation of zero-crossing counting on
%   this kind of data, not a bug to be tuned away — DO NOT trust this
%   function's output blindly. Always sanity-check n_cycles and
%   mean_freq_Hz against a visual beat count and/or the FFT-based
%   head_TBF/tail_TBF for the same trial. Passing ref_freq_Hz (e.g.
%   kine.tail_TBF) makes this function do that comparison automatically
%   and warn if the two disagree by more than 25%.
%
%   INPUTS
%     y             - [nFrames x 1] lateral (or DV) displacement time
%                     series. Should already be NaN-filled (see fill_nan
%                     in compute_kinematics.m) — this function does not
%                     itself handle missing data.
%     fps           - frames per second
%     min_freq      - (optional) Hz. Cycles implying a frequency below
%                     this are dropped. Default 0 (no floor).
%     max_freq      - (optional) Hz. Cycles implying a frequency above
%                     this are dropped. Default Inf (no ceiling).
%     ref_freq_Hz   - (optional) an independently-computed reference
%                     frequency (e.g. kine.tail_TBF) to cross-check
%                     against. Triggers a warning if mean_freq_Hz differs
%                     from it by more than 25%.
%
%   OUTPUT  cyc — struct:
%     .freqs_Hz            [nCycles x 1] frequency of each detected cycle
%     .periods_s           [nCycles x 1] period of each detected cycle
%     .cycle_start_frame   [nCycles x 1] fractional (interpolated) frame
%                          index where each cycle begins
%     .cycle_end_frame     [nCycles x 1] same, where each cycle ends
%     .n_cycles            scalar count
%     .mean_freq_Hz / .std_freq_Hz   NaN if n_cycles == 0
%     .detrend_window_frames         the moving-average window actually used

    if nargin < 3 || isempty(min_freq), min_freq = 0; end
    if nargin < 4 || isempty(max_freq), max_freq = Inf; end
    if nargin < 5, ref_freq_Hz = NaN; end

    y = y(:);
    n = numel(y);

    if all(isnan(y)) || n < 10
        cyc = empty_cycle_result();
        return;
    end

    % ---- Detrend: subtract a moving-average baseline before counting
    % zero-crossings, so slow drift isn't mistaken for the DC level of a
    % clean oscillation. Window = n/5 by default — a heuristic, not a
    % universally correct value; if your cycle counts look wrong, this is
    % the first thing to try adjusting (edit the constant below or copy
    % this function and expose it as a parameter for your use case). ----
    detrend_window = max(5, round(n/5));
    pad = floor(detrend_window/2);
    y_padded = [repmat(y(1),pad,1); y; repmat(y(end),pad,1)];
    kernel = ones(detrend_window,1) / detrend_window;
    trend_full = conv(fillmissing(y_padded,'linear'), kernel, 'valid');
    trend = trend_full(1:n);
    y_detrend = y - trend;

    % ---- Positive-going zero-crossings, sub-frame interpolated ----
    crossings = [];
    for k = 1:n-1
        yk = y_detrend(k); yk1 = y_detrend(k+1);
        if isnan(yk) || isnan(yk1), continue; end
        if yk <= 0 && yk1 > 0
            frac = -yk / (yk1 - yk);
            crossings(end+1) = k + frac; %#ok<AGROW>
        end
    end

    if numel(crossings) < 2
        cyc = empty_cycle_result();
        cyc.detrend_window_frames = detrend_window;
        return;
    end

    periods_frames = diff(crossings);
    periods_s      = periods_frames / fps;
    freqs_Hz       = 1 ./ periods_s;

    valid = freqs_Hz >= min_freq & freqs_Hz <= max_freq;

    starts = crossings(1:end-1)';
    ends   = crossings(2:end)';

    cyc.freqs_Hz          = freqs_Hz(valid)';
    cyc.periods_s         = periods_s(valid)';
    cyc.cycle_start_frame = starts(valid);
    cyc.cycle_end_frame   = ends(valid);
    cyc.n_cycles          = numel(cyc.freqs_Hz);
    cyc.detrend_window_frames = detrend_window;

    if cyc.n_cycles > 0
        cyc.mean_freq_Hz = mean(cyc.freqs_Hz);
        cyc.std_freq_Hz  = std(cyc.freqs_Hz);
    else
        cyc.mean_freq_Hz = NaN;
        cyc.std_freq_Hz  = NaN;
    end

    if ~isnan(ref_freq_Hz) && ~isnan(cyc.mean_freq_Hz) && ref_freq_Hz > 0
        rel_diff = abs(cyc.mean_freq_Hz - ref_freq_Hz) / ref_freq_Hz;
        if rel_diff > 0.25
            warning(['compute_cycle_frequencies: per-cycle mean (%.3f Hz, %d cycles) ' ...
                     'differs from the reference frequency (%.3f Hz) by %.0f%% — this is ' ...
                     'exactly the kind of disagreement flagged in this function''s docstring. ' ...
                     'Do a visual beat count on this trial before trusting either number.'], ...
                     cyc.mean_freq_Hz, cyc.n_cycles, ref_freq_Hz, 100*rel_diff);
        end
    end
end

function cyc = empty_cycle_result()
    cyc.freqs_Hz = []; cyc.periods_s = [];
    cyc.cycle_start_frame = []; cyc.cycle_end_frame = [];
    cyc.n_cycles = 0; cyc.mean_freq_Hz = NaN; cyc.std_freq_Hz = NaN;
    cyc.detrend_window_frames = NaN;
end
