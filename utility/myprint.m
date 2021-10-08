function str = myprint(v, precision)
    % convert v into a vector of numbers with 
    
    fmt0 = append('%.',num2str(precision),'g');
    if isempty(v)
        str = '\n';
    elseif isscalar(v)
        fmt = append(fmt0, '\n');
        str = sprintf(fmt, v);
    else
        N = numel(v);
        fmt = append(fmt0,' ');
        fmt = append(repmat(fmt, 1, N-1), fmt0, '\n');
        str = sprintf(fmt, v(:));
    end
end