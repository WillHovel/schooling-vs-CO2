function check_all_close(actual, expected, tol, name)
% CHECK_ALL_CLOSE  Elementwise comparison; reports the worst deviation.
%   NaN patterns must match exactly (NaN is the toolkit's "no data" value).
    actual = actual(:);
    expected = expected(:);
    if numel(actual) ~= numel(expected)
        error('check:failed', '%s: size mismatch (%d vs %d elements).', ...
              name, numel(actual), numel(expected));
    end
    nan_a = isnan(actual);
    nan_e = isnan(expected);
    if any(nan_a ~= nan_e)
        error('check:failed', '%s: NaN pattern mismatch (%d vs %d NaNs).', ...
              name, sum(nan_a), sum(nan_e));
    end
    both = ~nan_a;
    if ~any(both), return; end
    scale = max(1, abs(expected(both)));
    dev = abs(actual(both) - expected(both)) ./ scale;
    [mx, mi] = max(dev);
    if mx > tol
        error('check:failed', '%s: worst dev %.3g at element %d (got %.8g, expected %.8g, tol %.3g).', ...
              name, mx, mi, actual(both(mi)), expected(both(mi)), tol);
    end
end
