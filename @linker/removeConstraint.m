function removeConstraint(obj, name)
    % removeConstraint is a public method for the class LINKER. 
    % 
    %   removeConstraint(obj, name) remove the enforced constraint over the
    %   input variable.
    %       Input Arguments:
    %                   obj: A column linker object
    %                  name: Input of type char; semantic variable name.
    %                        See getParamsName(obj) for a list of valid
    %                        semantic names.
    %
    
    validateattributes(obj, {'linker'}, {'column'});
    validateattributes(name, {'char'}, {});
    assert(all([obj.varNum]==size(obj,1)), ...
        'MATLAB:linker:updateCompOrder',...
        "'varNum' must be the same with the length of input object.");
    
    varNames = {obj.name};
    idx = strcmp(name, varNames);
    assert(sum(idx)==1, 'MATLAB:linker:removeConstraint', ...
        "Unknown name '%s'. Valid variable names include:\n%s",...
        name, strjoin(varNames, ', '));
    
    % check if variable is dependent
    if strcmp(obj(idx).type, 'independent')
        warning('MATLAB:linker:removeConstrain',...
            '%s is already an independent variable.', name);
        return;
    end
    
    % remove constraint
    obj(idx).type = 'independent';
    obj(idx).constraint = '';
    obj(idx).dependence = [];    
    if obj(idx).bounded
        obj(idx).fval = @(x,~, u, l) (u-l)*cos(x)^2+l;
        obj(idx).fvalinv = @(p, u, l) acos(sqrt((p-l)/(u-l)));
        obj(idx).fgrad = @(x, u, l) -sin(2*x)*(u - l);
    else
        obj(idx).fval = @(x, ~, ~, ~) x;
        obj(idx).fvalinv = @(p, ~, ~) p;
        obj(idx).fgrad = @(~, ~, ~) 1;
    end
    obj(idx).ubound = @(~, u, ~) u;
    obj(idx).uboundgrad = @(~, ~, N) zeros(1,N);
    obj(idx).lbound = @(~, ~, l) l;
    obj(idx).lboundgrad = @(~, ~, N) zeros(1,N);
    obj.updateCompOrder();
end