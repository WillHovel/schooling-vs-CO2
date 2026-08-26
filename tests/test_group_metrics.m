function test_group_metrics()
% TEST_GROUP_METRICS  Ground-truth checks for the group-level metric
%                     functions: compute_polarization, compute_angle_to_flow,
%                     compute_distance_between_individuals.
%
%   Three hand-built fish with exactly-known snout/peduncle positions:
%     A: snout (1,0), ped (0,0)      -> heading 0 deg
%     B: snout (4,4), ped (3,4)      -> heading 0 deg, 5 units from A
%     C: snout (9,0), ped (9,-1)     -> heading 90 deg, 8 from A,
%                                       sqrt(41) from B
%   1. polarization of aligned school = 1 exactly; [0,90] -> sqrt(2)/2;
%      [0,0,90] -> sqrt(5)/3; frame with 1 fish -> NaN.
%   2. angle_to_flow: heading 90 with flow axis 0 -> 90; flow axis 90 ->
%      0; the +/-179 wrap case -> circular mean 180 (not 0), wrap-safe
%      range 2, circular std sqrt(-2*ln(cosd(1))) in deg.
%   3. distance: exact pairwise distances, symmetric matrix, NaN diagonal,
%      mean-NN and mean-pairwise formulas, cm_per_unit scaling.

    nF = 5;
    mk = @(name, sn, pd) struct('name', name, 'frames', (1:nF)', ...
        'point_names', {{'snout','ped'}}, ...
        'points', cat(2, repmat(reshape(sn, 1, 1, 2), nF, 1), ...
                          repmat(reshape(pd, 1, 1, 2), nF, 1)), ...
        'has_z', false);
    fpA = mk('A', [1 0],   [0 0]);
    fpB = mk('B', [4 4],   [3 4]);
    fpC = mk('C', [9 0],   [9 -1]);

    % ---- 1) Polarization ----
    pol2 = compute_polarization([fpA fpB], 'snout', 'ped');
    check_all_close(pol2.heading_deg, zeros(nF,2), 1e-9, 'both fish heading 0');
    check_all_close(pol2.polarization, ones(nF,1), 1e-12, 'aligned school -> pol = 1');
    check_close(pol2.mean_polarization, 1, 1e-12, 'mean polarization = 1');
    check_all_close(pol2.n_fish_present, 2*ones(nF,1), 0, 'both present every frame');

    polAC = compute_polarization([fpA fpC], 'snout', 'ped');
    check_close(polAC.polarization(1), sqrt(2)/2, 1e-12, '[0,90] -> pol = sqrt(2)/2');
    check_close(polAC.heading_deg(1,2), 90, 1e-9, 'fish C heading = 90');

    pol3 = compute_polarization([fpA fpB fpC], 'snout', 'ped');
    check_close(pol3.polarization(1), sqrt(5)/3, 1e-12, '[0,0,90] -> pol = sqrt(5)/3');

    fpB2 = fpB; fpB2.points(3, 1, :) = NaN;
    polN = compute_polarization([fpA fpB2], 'snout', 'ped');
    check(isnan(polN.polarization(3)), 'frame with 1 fish -> pol NaN');
    check(polN.n_fish_present(3) == 1, 'n_present = 1', 'got %d', polN.n_fish_present(3));
    check_close(polN.mean_polarization, 1, 1e-12, 'mean over valid frames = 1');

    % ---- 2) Angle to flow ----
    angC = compute_angle_to_flow(fpC, 'snout', 'ped', 0);
    check_all_close(angC.angle_to_flow_deg, 90*ones(nF,1), 1e-9, 'heading 90, flow 0 -> 90');
    check_close(angC.mean_angle_to_flow_deg, 90, 1e-9, 'circular mean = 90');
    check_close(angC.std_angle_to_flow_deg, 0, 1e-9, 'circular std = 0');
    check_close(angC.range_angle_to_flow_deg, 0, 1e-9, 'range = 0');
    check(angC.n_valid_frames == nF, 'n_valid = 5', 'got %d', angC.n_valid_frames);

    angC90 = compute_angle_to_flow(fpC, 'snout', 'ped', 90);
    check_all_close(angC90.angle_to_flow_deg, zeros(nF,1), 1e-9, 'flow axis 90 -> angle 0');

    fpW = mk('W', [cosd(179) sind(179)], [0 0]);
    fpW.points(2, 1, :) = [cosd(-179) sind(-179)];
    % keep only the two wrap frames: with equal +179/-179 counts the
    % circular mean is exactly 180 (5 frames with a 4:1 split would give
    % 179.4, which is a count artifact, not a wrap bug)
    fpW.points = fpW.points(1:2, :, :);
    fpW.frames = fpW.frames(1:2);
    angW = compute_angle_to_flow(fpW, 'snout', 'ped');
    check_close(angW.heading_deg(1), 179, 1e-9, 'wrap frame 1 = 179');
    check_close(angW.heading_deg(2), -179, 1e-9, 'wrap frame 2 = -179');
    check_close(angW.mean_angle_to_flow_deg, 180, 1e-9, 'circular mean = 180, not 0');
    check_close(angW.std_angle_to_flow_deg, sqrt(-2*log(cosd(1)))*180/pi, 1e-9, ...
                'circular std = sqrt(-2 ln R) in deg');
    check_close(angW.range_angle_to_flow_deg, 2, 1e-9, 'wrap-safe range = 2');

    % ---- 3) Distance between individuals ----
    dAB = 5; dAC = 8; dBC = sqrt(41);
    dist2 = compute_distance_between_individuals([fpA fpB fpC], 'snout', 2);
    check_close(dist2.cm_per_unit, 2, 0, 'cm_per_unit stored');
    check_close(dist2.pairwise_dist(1,1,2), 2*dAB, 1e-9, 'd(A,B) = 2*5');
    check_close(dist2.pairwise_dist(1,2,1), 2*dAB, 1e-9, 'matrix symmetric');
    check_close(dist2.pairwise_dist(1,1,3), 2*dAC, 1e-9, 'd(A,C) = 2*8');
    check_close(dist2.pairwise_dist(1,2,3), 2*dBC, 1e-9, 'd(B,C) = 2*sqrt(41)');
    check(isnan(dist2.pairwise_dist(1,1,1)), 'diagonal NaN');
    check_close(dist2.mean_nn_dist(1), 2*(dAB + dAB + dBC)/3, 1e-9, ...
                'mean NN = (10+10+2sqrt(41))/3');
    check_close(dist2.mean_pairwise_dist(1), 2*(dAB + dAC + dBC)/3, 1e-9, ...
                'mean pairwise = (10+16+2sqrt(41))/3');
    check_close(dist2.overall_mean_nn_dist, 2*(dAB + dAB + dBC)/3, 1e-9, 'overall NN');
    check(isequal(dist2.fish_names, {'A','B','C'}), 'fish_names order kept');

    dist1 = compute_distance_between_individuals([fpA fpB], 'snout');
    check_close(dist1.pairwise_dist(1,1,2), dAB, 1e-9, 'default cm_per_unit = 1');
    check_close(dist1.overall_mean_nn_dist, dAB, 1e-9, 'two-fish NN = dAB');
end
