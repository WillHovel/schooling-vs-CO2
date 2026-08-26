function test_cycle_frequencies()
% TEST_CYCLE_FREQUENCIES  Ground-truth checks for compute_cycle_frequencies.m.
%
%   A pure 2 Hz sinusoid sampled at 100 fps for 4 s: the detrended signal
%   is still an exact 2 Hz sinusoid in the interior (a moving average of a
%   sinusoid is a sinusoid), so per-cycle frequencies must measure 2 Hz.
%   Interior cycles are exact; the first and last cycles touch the padded
%   edges of the detrend window, so the mean is checked to 3% and the
%   spread to 0.05 Hz. 7 full cycles fit in 4 s.
%   Then the min_freq/max_freq gates, the ref_freq cross-check warning,
%   and the degenerate inputs (all-NaN, too short, constant).

    fps = 100; nF = 400;
    t = (0:nF-1)'/fps;
    y = sin(2*pi*2*t);

    cyc = compute_cycle_frequencies(y, fps);
    check(cyc.n_cycles == 7, '2 Hz over 4 s -> 7 full cycles', 'got %d', cyc.n_cycles);
    check(numel(cyc.freqs_Hz) == 7, 'per-cycle freqs exported');
    check(numel(cyc.periods_s) == 7, 'per-cycle periods exported');
    check_close(cyc.mean_freq_Hz, 2, 3e-2, 'mean freq = 2 Hz');
    check_close(cyc.std_freq_Hz, 0, 5e-2, 'clean sinusoid: per-cycle spread ~0');
    check_close(cyc.periods_s(4), 0.5, 5e-3, 'interior period = 0.5 s');

    % ---- min/max frequency gates ----
    cycHi = compute_cycle_frequencies(y, fps, 2.5);
    check(cycHi.n_cycles == 0, 'min_freq 2.5 gate drops all cycles', 'got %d', cycHi.n_cycles);
    check(isnan(cycHi.mean_freq_Hz), 'gated-out mean is NaN');
    cycLo = compute_cycle_frequencies(y, fps, 0, 1);
    check(cycLo.n_cycles == 0, 'max_freq 1 gate drops all cycles', 'got %d', cycLo.n_cycles);

    % ---- ref_freq cross-check warning ----
    lastwarn('');
    cycW = compute_cycle_frequencies(y, fps, 0, Inf, 1);
    [wmsg, ~] = lastwarn;
    check(contains(wmsg, 'differs from the reference'), 'ref_freq mismatch warns', ...
          'lastwarn = "%s"', wmsg);
    check(numel(cycW.periods_s) == 7, 'ref_freq warning leaves cycles unchanged');

    % ---- degenerate inputs ----
    cycN = compute_cycle_frequencies(NaN(400,1), fps);
    check(cycN.n_cycles == 0 && isnan(cycN.mean_freq_Hz), 'all-NaN -> empty result');
    cycS = compute_cycle_frequencies(zeros(7,1), fps);
    check(cycS.n_cycles == 0, 'too-short input -> empty result');
    cycC = compute_cycle_frequencies(ones(400,1), fps);
    check(cycC.n_cycles == 0, 'constant signal -> no crossings');
end
