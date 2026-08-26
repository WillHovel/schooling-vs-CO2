function test_app_build()
% TEST_APP_BUILD  Smoke test: the Kinemetrix GUI must construct and
% delete cleanly. FishKinematicsApp has no output argument (plain
% function that builds a uifigure named Kinemetrix), so the figure is
% looked up by name.

    FishKinematicsApp();
    h = findall(groot, 'Type', 'figure', 'Name', 'Kinemetrix');
    check(~isempty(h), 'Kinemetrix figure created on startup');
    delete(h);
end
