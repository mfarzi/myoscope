function n = getVarNum(obj, filter)
    % return the number of variables
    if nargin==1
        filter = 'all';
    end
    
    switch filter
        case 'all'
            n = size(obj, 1);
        case 'dummy' % a copy of another free parameter
            isDummy = strcmp({obj.type}', 'dummy');
            n = sum(isDummy);
        case 'constant' % constant value
            isConstant = strcmp({obj.type}', 'constant');
            n = sum(isConstant);
        case 'free' % optimisable variables
            isDummy = strcmp({obj.type}', 'dummy');
            isConstant = strcmp({obj.type}', 'constant');
            n = sum(not(isDummy|isConstant));
        case 'independent' %parameters which are not dependent on others
            isIndependent = strcmp({obj.type}', 'independent');
            n = sum(isIndependent);
        case 'dependent' % parameters which depend on other variables
            isDependent = strcmp({obj.type}', 'dependent');
            n = sum(isDependent);
        otherwise
            error('MATLAB:linker:getVarNum','Unknown filter.');
    end
end
            