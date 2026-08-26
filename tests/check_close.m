function check_close(actual, expected, tol, name)
% CHECK_CLOSE  Compare a scalar metric to its expected value.
%   tol is absolute for |expected| <= 1, otherwise relative to |expected|.
%   NaN expectations require NaN actuals (the toolkit's "no data" convention).
    if isnan(expected)
        if ~isnan(actual)
            error('check:failed', '%s: expected NaN, got %.6g.', name, actual);
        end
        return;
    end
    if isnan(actual)
        error('check:failed', '%s: got NaN, expected %.6g +/- %.3g.', name, expected, tol);
    end
    scale = max(1, abs(expected));
    dev = abs(actual - expected) / scale;
    if dev > tol
        error('check:failed', '%s: got %.8g, expected %.8g (dev %.3g > tol %.3g).', ...
              name, actual, expected, dev, tol);
    end
end
