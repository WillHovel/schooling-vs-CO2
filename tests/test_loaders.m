function test_loaders()
% TEST_LOADERS  Ground-truth checks for the five CSV loaders.
%
%   Writes small synthetic CSVs with write_csv and checks that every
%   loader reproduces the exact values written, the right schema
%   (point_names / has_z / format / frames), anatomical labels where
%   applicable, and the special paths (compound points, all-missing text
%   columns, dual-camera pectoral phase, p_cutoff likelihood filtering,
%   CURVES pre-transformed passthrough).

    d = tempname;
    mkdir(d);
    cleanup = onCleanup(@() rmdir(d, 's'));
    f = @(nm) fullfile(d, nm);

    % ==================== FORMAT A: DLC-style Fish1_P1_x ====================
    rowsA = [ {'frame','Fish1_P1_x','Fish1_P1_y','Fish1_P2_x','Fish1_P2_y'}
              num2cell([(1:5)', (1:5)'*1, (1:5)'*2, (1:5)'*3, (1:5)'*4]) ];
    write_csv(f('fmtA.csv'), rowsA);
    fpA = load_fish_points(f('fmtA.csv'));
    check(strcmp(fpA.format, 'DLC'), 'A: format = DLC');
    check(~fpA.has_z, 'A: has_z = false');
    check(isequal(fpA.frames, (1:5)'), 'A: frames = frame column');
    check(isequal(fpA.point_names, {'P1','P2'}), 'A: point_names = P1,P2');
    check_all_close(fpA.points(:, :, 1), [(1:5)', (1:5)'*3], 0, 'A: X values exact');
    check_all_close(fpA.points(:, :, 2), [(1:5)'*2, (1:5)'*4], 0, 'A: Y values exact');

    % Format A 3D + multi-fish
    rowsA3 = [ {'frame','Fish1_P1_x','Fish1_P1_y','Fish1_P1_z','Fish2_P1_x','Fish2_P1_y','Fish2_P1_z'}
               num2cell([(1:4)', (1:4)', (1:4)'*2, (1:4)'*3, (1:4)'*10, (1:4)'*20, (1:4)'*30]) ];
    write_csv(f('fmtA3.csv'), rowsA3);
    fpA3 = load_fish_points(f('fmtA3.csv'));
    check(numel(fpA3) == 2, 'A3: two fish found', 'got %d', numel(fpA3));
    check(strcmp(fpA3(2).name, 'Fish2'), 'A3: second fish name = Fish2');
    check(fpA3(1).has_z, 'A3: has_z = true');
    check_all_close(fpA3(1).points(:, 1, 3), (1:4)'*3, 0, 'A3: Z values exact');
    check_all_close(fpA3(2).points(:, 1, 1), (1:4)'*10, 0, 'A3: Fish2 X exact');

    % ==================== FORMAT C: numbered ptN_X ====================
    rowsC = [ {'pt13_X','pt13_Y','pt1_X','pt1_Y','pt1_Z','pt2_X','pt2_Y','pt2_Z'}
              num2cell([(1:4)', (1:4)'*9, (1:4)'*2, (1:4)'*3, (1:4)'*4, ...
                        (1:4)'*5, (1:4)'*6, (1:4)'*7]) ];
    write_csv(f('fmtC.csv'), rowsC);
    fpC = load_fish_points(f('fmtC.csv'));
    check(strcmp(fpC.format, 'numbered'), 'C: format = numbered');
    check(fpC.has_z, 'C: has_z = true');
    check(isequal(fpC.point_names, {'pt1_pect_base_R','pt2_pect_tip_R','pt13'}), ...
          'C: anatomical labels + unknown point stays raw', 'got {%s}', strjoin(fpC.point_names, ','));
    check_all_close(fpC.points(:, 1, 1), (1:4)'*2, 0, 'C: pt1 X exact (sorted numerically)');
    check_all_close(fpC.points(:, 3, 2), (1:4)'*9, 0, 'C: pt13 Y exact');

    % Format C 2D
    rowsC2 = [ {'pt1_X','pt1_Y','pt2_X','pt2_Y'}
               num2cell([(1:4)', (1:4)', (1:4)'*2, (1:4)'*3]) ];
    write_csv(f('fmtC2.csv'), rowsC2);
    fpC2 = load_fish_points(f('fmtC2.csv'));
    check(~fpC2.has_z, 'C2: has_z = false');

    % ==================== FORMAT B: named columns ====================
    rowsB = [ {'snout_X','snout_Y','snout_Z','tail_X','tail_Y','tail_Z'}
              num2cell([(1:4)', (1:4)'*2, (1:4)'*3, (1:4)'*4, (1:4)'*5, (1:4)'*6]) ];
    write_csv(f('fmtB.csv'), rowsB);
    fpB = load_fish_points_named(f('fmtB.csv'), {'snout','tail'});
    check(fpB.has_z, 'B: has_z = true');
    check(isequal(fpB.point_names, {'snout','tail'}), 'B: selection order kept');
    check_all_close(fpB.points(:, 1, 1), (1:4)', 0, 'B: snout X exact');
    check_all_close(fpB.points(:, 2, 3), (1:4)'*6, 0, 'B: tail Z exact');

    % point_order permutation
    fpBp = load_fish_points_named(f('fmtB.csv'), {'snout','tail'}, [2 1]);
    check(isequal(fpBp.point_names, {'tail','snout'}), 'B: point_order permutation');
    check_all_close(fpBp.points(:, 1, 1), (1:4)'*4, 0, 'B: permuted tail X exact');

    % compound point (mean of a pair)
    rowsB2 = [ {'snout_L_X','snout_L_Y','snout_R_X','snout_R_Y','tail_X','tail_Y'}
               num2cell([(1:4)', (1:4)'*2, (1:4)'*3, (1:4)'*4, (1:4)'*5, (1:4)'*6]) ];
    write_csv(f('fmtB2.csv'), rowsB2);
    fpB2 = load_fish_points_named(f('fmtB2.csv'), {{'snout_L','snout_R'}, 'tail'});
    check(isequal(fpB2.point_names, {'snout_L+snout_R','tail'}), 'B: compound label');
    check_all_close(fpB2.points(:, 1, 1), (1:4)'*2, 0, 'B: compound X = mean of pair');

    % no-selection mode lists available names, returns empty points
    fpB0 = load_fish_points_named(f('fmtB.csv'));
    check(isempty(fpB0.points), 'B: no-selection returns points = []');
    check(isequal(fpB0.point_names, {'snout','tail'}), 'B: no-selection lists names');

    % all-'NA' column: must come out as clean NaN (imported as text or
    % numeric depending on readtable's missing-value handling — either way
    % the column must not crash and must be all-NaN)
    rowsB3 = [ {'snout_X','snout_Y','occl_X','occl_Y'}
               num2cell([(1:4)', (1:4)'*2, nan(4,1), nan(4,1)]) ];
    rowsB3{2,3} = 'NA'; rowsB3{3,3} = 'NA'; rowsB3{4,3} = 'NA'; rowsB3{5,3} = 'NA';
    rowsB3{2,4} = 'NA'; rowsB3{3,4} = 'NA'; rowsB3{4,4} = 'NA'; rowsB3{5,4} = 'NA';
    write_csv(f('fmtB3.csv'), rowsB3);
    fpB3 = load_fish_points_named(f('fmtB3.csv'), {'snout','occl'});
    check(all(isnan(fpB3.points(:, 2, :)), 'all'), 'B: all-NA column -> all-NaN');
    check_all_close(fpB3.points(:, 1, 1), (1:4)', 0, 'B: good column beside all-NA one intact');

    % guaranteed-text column (never a missing indicator): warning path
    rowsB4 = rowsB3;
    rowsB4{2,3} = 'missed'; rowsB4{3,3} = 'missed'; rowsB4{4,3} = 'missed'; rowsB4{5,3} = 'missed';
    write_csv(f('fmtB4.csv'), rowsB4);
    lastwarn('');
    fpB4 = load_fish_points_named(f('fmtB4.csv'), {'snout','occl'});
    [wmsg, ~] = lastwarn;
    check(all(isnan(fpB4.points(:, 2, 1))), 'B: text column -> all-NaN');
    check(contains(wmsg, '100%'), 'B: text column warns about 100%% missing', ...
          'lastwarn = "%s"', wmsg);

    % ==================== FORMAT D: dual camera ====================
    nF = 100; t = (0:nF-1)';
    y2 = sin(2*pi*2*t/nF);
    hdrD = {'pt1_cam1_X','pt1_cam1_Y','pt1_cam2_X','pt1_cam2_Y', ...
            'pt2_cam1_X','pt2_cam1_Y','pt12_cam1_X','pt12_cam1_Y'};
    rowsD = [hdrD; num2cell([(1:nF)', (1:nF)'*0.1, (1:nF)'*0.2, (1:nF)'*0.3, ...
                             (1:nF)'*0.4, y2, (1:nF)'*0.5, y2])];
    write_csv(f('fmtD_in.csv'), rowsD);
    fpD = load_fish_points_named(f('fmtD_in.csv'));
    check(strcmp(fpD.format, 'dual_camera'), 'D: format = dual_camera');
    check(isequal(fpD.cam_names, {'cam1','cam2'}), 'D: cam_names sorted');
    check(isequal(fpD.point_names, {'pt1_pect_base_R','pt2_pect_tip_R','pt12_pect_base_L'}), ...
          'D: anatomical labels', 'got {%s}', strjoin(fpD.point_names, ','));
    check_all_close(fpD.points(:, 2, 2), y2, 1e-9, 'D: points = first camera values');
    check(strcmp(fpD.pect_phase_result.classification, 'In-phase'), 'D: in-phase classification');
    check_close(fpD.pect_phase_result.phase_shift_deg, 0, 1e-6, 'D: in-phase shift = 0');
    check(fpD.pect_phase_result.n_valid == nF, 'D: n_valid = nFrames');

    rowsD2 = rowsD;
    rowsD2(2:end, 8) = num2cell(-y2);
    write_csv(f('fmtD_anti.csv'), rowsD2);
    fpD2 = load_fish_points_named(f('fmtD_anti.csv'));
    check(strcmp(fpD2.pect_phase_result.classification, 'Antiphase'), 'D: antiphase classification');
    check_close(fpD2.pect_phase_result.phase_shift_deg, 180, 1e-6, 'D: antiphase shift = 180');

    % ==================== FORMAT E: CURVES ====================
    nSt = 19; mm = 60:10:240;
    hdrE = cell(1, 2*nSt);
    hdrE(1:2:end) = num2cell(mm);                 % odd columns: station mm
    hdrE(2:2:end) = {NaN};                        % even columns: NaN
    Xv = repmat(linspace(0, 1, nSt), 4, 1);
    Yv = 0.05 * repmat(sin(linspace(0, 2*pi, nSt)), 4, 1);
    dataE = NaN(4, 2*nSt);
    dataE(:, 1:2:end) = Xv;
    dataE(:, 2:2:end) = Yv;
    rowsE = [hdrE; num2cell(dataE); num2cell(nan(1, 2*nSt))];   % trailing blank row
    write_csv(f('fmtE.csv'), rowsE);
    fpE = load_fish_curves(f('fmtE.csv'));
    check(fpE.pre_transformed, 'E: pre_transformed = true');
    check_all_close(fpE.body_pos_mm, mm, 0, 'E: station positions 60:10:240');
    check(isequal(fpE.point_names{1}, 's60mm') && isequal(fpE.point_names{end}, 's240mm'), ...
          'E: station labels');
    check(size(fpE.points, 1) == 4, 'E: trailing all-NaN row dropped', ...
          'got %d', size(fpE.points, 1));
    check_all_close(fpE.X, Xv, 1e-9, 'E: X = normalized station positions');
    check_all_close(fpE.Y, Yv, 1e-9, 'E: Y exact');

    % ==================== DLC multi-animal ====================
    rowsM = { 'scorer','n1','n1','n1','n1','n1','n1','n1','n1','n1','n1','n1','n1'
              'individuals','indivA','indivA','indivA','indivA','indivA','indivA', ...
                            'indivB','indivB','indivB','indivB','indivB','indivB'
              'bodyparts','snout','snout','snout','mid','mid','mid', ...
                          'snout','snout','snout','mid','mid','mid'
              'coords','x','y','likelihood','x','y','likelihood', ...
                       'x','y','likelihood','x','y','likelihood'
              1, 1, 2, 0.9, 3, 4, 0.9, 10, 20, 0.9, 30, 40, 0.9
              2, 5, 6, 0.5, 7, 8, 0.9, 50, 60, 0.9, 70, 80, 0.9
              3, 9, 10, 0.9, 11, 12, 0.9, 90, 100, 0.9, 110, 120, 0.9 };
    write_csv(f('fmtM.csv'), rowsM);
    fpM = load_fish_points_dlc_multianimal(f('fmtM.csv'), 0.6);
    check(numel(fpM) == 2, 'M: two individuals', 'got %d', numel(fpM));
    check(strcmp(fpM(1).name, 'indivA') && strcmp(fpM(2).name, 'indivB'), 'M: names');
    check(isequal(fpM(1).point_names, {'snout','mid'}), 'M: bodypart order');
    check_all_close(fpM(1).points(:, 1, 1), [1; NaN; 9], 0, 'M: low-likelihood frame NaN (p_cutoff)');
    check_all_close(fpM(1).points(:, 2, 2), [4; 8; 12], 0, 'M: second bodypart exact');
    check_all_close(fpM(2).points(:, 1, 1), [10; 50; 90], 0, 'M: indivB untouched by cutoff');
    check_all_close(fpM(1).frames, [1; 2; 3], 0, 'M: frames from column 1');
end
