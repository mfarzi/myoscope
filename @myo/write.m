function write(fileId, fieldname, v)
    % write is public method for class myo
    nrows = size(v, 1);
    
    if isnumeric(v)
        fprintf(fileId, '#$%s %d numeric\n', fieldname, nrows);
        for i=1:rows
            myo.print(fileId, v(i,:), '%0.10g');
        end
    else
        fprintf(fileId, '#$%s %d string\n', fieldname, nrows);
    end
    
end