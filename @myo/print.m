function print(fileId, v, fmt0)
    % convert v into a vector of numbers with 
    
    %fmt0 = append('%.',num2str(precision),'g');
    if isempty(v)
        fprintf(fileId, '%s', '\n');
    elseif isscalar(v)
        fmt = append(fmt0, '\n');
        fprintf(fileId, fmt, v);
    else
        N = numel(v);
        fmt = append(fmt0,' ');
        fmt = append(repmat(fmt, 1, N-1), fmt0, '\n');
        fprintf(fileId, fmt, v(:));
    end
end