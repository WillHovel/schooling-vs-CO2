function check(cond, name, varargin)
% CHECK  Assertion helper for the Kinemetrix test suite.
%   check(cond, name) errors with name if cond is false.
%   check(cond, name, fmt, ...) errors with name plus sprintf(fmt, ...).
% Throws (caught and reported by run_all_tests.m) rather than using
% assert() so failures carry the metric name and values.
    if ~isscalar(cond) || ~islogical(cond)
        error('check:condition', 'check(%s): condition must be a scalar logical.', name);
    end
    if ~cond
        if isempty(varargin)
            error('check:failed', '%s: check failed.', name);
        else
            error('check:failed', '%s: check failed. %s', name, sprintf(varargin{:}));
        end
    end
end
