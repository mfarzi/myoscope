function P = map(obj, X)  
    % map is a public method for the class linker
    %
    %   map(obj, X) maps unconstrained optimisation variables X into model 
    %   parameters P
    %
%     validateattributes(obj, {'linker'}, {'column'});
%     validateattributes(X, {'numeric'}, {'size', [obj.getVarNum('free'),1]},...
%         'MATLAB:linker:map', 'x');
    assert(not(any(isnan(X))), 'MATLAB:linker:map',...
       'Input has NaN values.');
    
    % initialise p with zeros
    N = size(obj, 1);
    P = zeros(N, 1);
    
    % add zero in place of dummy or constant variables
    X = zerofill(obj, X);
    
    % map scalar x to scalar p in order
    [~,index] = obj.getCompOrder();
    for n = index'
        u = obj(n).ubound(P(obj(n).dependence), obj(n).maximum, obj(n).minimum);
        l = obj(n).lbound(P(obj(n).dependence), obj(n).maximum, obj(n).minimum);
        P(n) = obj(n).fval(X(n), P(obj(n).dependence), u, l);
%         P(n) = obj(n).scalarMap(X(n), P);
    end
end