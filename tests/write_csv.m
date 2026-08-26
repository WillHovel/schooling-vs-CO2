function write_csv(fname, rows)
% WRITE_CSV  Write a cell matrix (strings / scalar numbers) as a CSV file.
%   rows: [nRows x nCols] cell array, all rows the same width.
    fid = fopen(fname, 'w');
    if fid < 0
        error('write_csv: cannot open %s for writing.', fname);
    end
    c = onCleanup(@() fclose(fid));
    for r = 1:size(rows, 1)
        for c_ = 1:size(rows, 2)
            v = rows{r, c_};
            if c_ > 1, fprintf(fid, ','); end
            if ischar(v) || isstring(v)
                fprintf(fid, '%s', char(v));
            elseif isscalar(v) && isnumeric(v)
                fprintf(fid, '%.10g', v);
            else
                error('write_csv: cell (%d,%d) must be a scalar or string.', r, c_);
            end
        end
        fprintf(fid, '\n');
    end
end
