function X = invmap(obj, P)
    % invmap is a public method for the class linker
    %
    %   invmap(obj, P) maps model parameters P to unconstrained 
    %   optimisation variables X
    %
%     validateattributes(obj, {'linker'}, {'column'});
%     validateattributes(P, {'numeric'}, {'size', [obj.getVarNum('all'),1]},...
%         'MATLAB:linker:invmap', 'p');
    assert(not(any(isnan(P))), 'MATLAB:linker:invmap',...
        'Input has NaN values.');

    
    % initialise X with zeros
    N = size(obj, 1);
    X = zeros(N, 1);

    % map scalar p to scalar x in order
    [~,index] = obj.getCompOrder();
    for n = index'
        u = obj(n).ubound(P(obj(n).dependence), obj(n).maximum, obj(n).minimum);
        l = obj(n).lbound(P(obj(n).dependence), obj(n).maximum, obj(n).minimum);
        X(n) = obj(n).fvalinv(P(n), u, l);
%         X(n) = obj(n).scalarInvmap(P(n), P);
    end

    % remove dummy zeros values
    isConstant = strcmp({obj.type}, 'constant');
    isDummy = strcmp({obj.type}, 'dummy');
    X(isConstant|isDummy) = [];
end