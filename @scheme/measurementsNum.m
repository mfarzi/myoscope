function n = measurementsNum(obj,str)
    % measurementsNum is a public method for class scheme
    %
    % return the number of measurements that satisfies the input constraint
    % 
    
    if nargin==1
        str = 'all';
    end
    
    switch str
        case 'all'
            n = size(obj.ghat, 1);
        case 'independent'
            n = 1;
    end
end