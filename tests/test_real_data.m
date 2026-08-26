function test_real_data()
% TEST_REAL_DATA  Regression checks against the manually-validated shark
% trials in Analysis/ and Shoval_Analysis/.
%
%   1. WSBS_Swim_1.5BL_Shark03xyzpts.csv: the known-good values visually
%      verified on this machine — head/tail/spline TBF = 1.8349 Hz,
%      wavelength in (1.1, 1.35) BL, wave_speed = 2.2296 BL/s.
%   2. WSBS_PSwalk_0.25BL_Shark01xyzpts.csv vs Shoval's manual tail-beat
%      peaks in Shoval_Analysis/TB_peaks.csv: the FFT-based tail_TBF must
%      stay within 15% of the mean manual peak-to-peak interval.
%
%   Skips (with a warning) if the data files are absent.

    here = fileparts(mfilename('fullpath'));
    root = fileparts(here);
    fps = 100;
    min_freq = 0.5;
    names = {'BP_1_SnoutML','BP_2_PectoralML','BP_3_PelvicML', ...
             'BP_4_AnalML','BP_5_CaudalML'};

    % ---- 1) Swim trial (known-good values) ----
    swim_csv = fullfile(root, 'Analysis', 'AllEdited_CSVs', 'AllEdited_CSVs', ...
                        'WSBS_Swim_1.5BL_Shark03xyzpts.csv');
    if ~isfile(swim_csv)
        warning('test_real_data: %s not found — skipping swim-trial checks.', swim_csv);
    else
        fp = load_fish_points_named(swim_csv, names);
        tf = transform_fish(fp);
        k = compute_kinematics(tf, fps, min_freq);
        check_close(k.head_TBF, 1.8349, 5e-2, 'swim trial head_TBF (known-good)');
        check_close(k.tail_TBF, 1.8349, 5e-2, 'swim trial tail_TBF');
        check_close(k.spline_freq_Hz, 1.8349, 5e-2, 'swim trial spline_freq');
        check(k.wavelength > 1.1 && k.wavelength < 1.35, ...
              'swim trial wavelength in (1.1, 1.35)', 'got %.4f', k.wavelength);
        check_close(k.wave_speed_BL_s, 2.2296, 1.5e-1, 'swim trial wave speed');
    end

    % ---- 2) Walking trial vs Shoval's manual peaks ----
    peaks_csv = fullfile(root, 'Shoval_Analysis', 'TB_peaks.csv');
    walk_csv  = fullfile(root, 'Analysis', 'AllEdited_CSVs', 'AllEdited_CSVs', ...
                         'WSBS_PSwalk_0.25BL_Shark01xyzpts.csv');
    if ~isfile(peaks_csv) || ~isfile(walk_csv)
        warning('test_real_data: Shoval cross-check files missing — skipping.');
        return;
    end
    T = readtable(peaks_csv, 'Delimiter', ',', 'VariableNamingRule', 'preserve');
    sel = startsWith(T.source_file, 'WSBS_PSwalk_0.25BL_Shark01xyzpts');
    peaks = T.Peak_Frames(sel);
    check(numel(peaks) >= 3, 'TB_peaks has Shark01 rows', 'got %d', numel(peaks));
    tbf_manual = fps / mean(diff(peaks));

    fpw = load_fish_points_named(walk_csv, names);
    tfw = transform_fish(fpw);
    kw = compute_kinematics(tfw, fps, min_freq);
    fprintf('  [manual TBF = %.4f Hz from %d peaks, computed tail_TBF = %.4f Hz]\n', ...
            tbf_manual, numel(peaks), kw.tail_TBF);
    check_close(kw.tail_TBF, tbf_manual, 1.5e-1, 'walking trial TBF vs manual peaks');
end
