function p = init(obj)
    dim = size(obj);
    if dim(2)>1 || length(dim) > 2
        error(['Linker object should be either ', ...
               'scalar or a column vector']);
    end

    p = zeros(dim);
    x = zeros(dim);
    for n = 1:dim(1)
        U = min(obj(n).upperBound, obj(n).initRange(2));
        L = max(obj(n).lowerBound, obj(n).initRange(1));
        p(n) = rand(1)*(U - L) + L;
        x(n) = obj(n).scalarInvLink(p(n));
    end

    % remove dummy or constant values
    isConstantOrDummy = obj.getConstant | obj.getDummy;
    x(isConstantOrDummy) = [];

    p = obj.link(x); 
end % of init