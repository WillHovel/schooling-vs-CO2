function run_all_tests()
% RUN_ALL_TESTS  Master runner for the Kinemetrix ground-truth test suite.
%
%   Adds the repo root and tests/ to the path, runs every test_*.m in
%   try/catch, prints PASS/FAIL per file plus the first error and its
%   stack, and errors out (non-zero exit) if any file failed — so it works
%   headlessly:
%       matlab -batch "cd('C:\Users\willh\Desktop\fish_analysis_v2'); run_all_tests"

    here = fileparts(mfilename('fullpath'));
    root = fileparts(here);
    addpath(root, here);

    tests = { @test_loaders, @test_transform_fish, @test_apply_body_transform, ...
              @test_filter_dlc_jumps, @test_compute_kinematics, ...
              @test_body_extended, @test_group_metrics, @test_cycle_frequencies, ...
              @test_fin_kinematics, @test_girdle_kinematics, @test_real_data, ...
              @test_app_build };

    nPass = 0; nFail = 0; failures = {};
    for i = 1:numel(tests)
        name = func2str(tests{i});
        fprintf('\n===== %s =====\n', name);
        try
            tests{i}();
            nPass = nPass + 1;
            fprintf('===== %s: PASS =====\n', name);
        catch err
            nFail = nFail + 1;
            failures{end+1} = name; %#ok<AGROW>
            fprintf('===== %s: FAIL =====\n', name);
            fprintf('  %s\n', err.message);
            for k = 1:numel(err.stack)
                fprintf('    at %s line %d\n', err.stack(k).name, err.stack(k).line);
            end
        end
    end

    fprintf('\n%d/%d test files passed.\n', nPass, nPass + nFail);
    if nFail > 0
        error('run_all_tests: %d test file(s) failed: %s', nFail, strjoin(failures, ', '));
    end
end
