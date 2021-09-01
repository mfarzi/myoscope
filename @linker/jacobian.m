function jac = jacobian(obj, x)
    % jacobian is a public method for the class linker
    %
    %   grad(obj, x) return the jacobian of forward mapping functions wrt
    %   optimisation variables x.
    
%     validateattributes(obj, {'linker'}, {'column'});
%     assert(all([obj.varNum]==size(obj,1)), ...
%         'MATLAB:linker:updateCompOrder',...
%         "'varNum' must be the same with the length of input object.");
    
    N = size(obj,1);
    p = obj.map(x);
    jac = eye(N);
    % add zero in place of dummy or constant variables
    x = zerofill(obj, x);
    

    [~,thisCompOrder] = obj.getCompOrder();
    for n = thisCompOrder'
        jac(n,:) = obj(n).grad(x, p, jac);
    end

    % remove dummy zeros values
    isDummy = strcmp({obj.type}, 'dummy');
    isConstant = strcmp({obj.type}, 'constant');
    jac(:,isDummy|isConstant) = [];
end % of jacobian

% function g = grad(obj, x, n, p, thisJac)
%     nObj = length(p);
% 
%     switch obj.type
%         case 'cos'
%             lBound = obj.lowerBound;
%             if obj.linked == 0
%                 uBound = obj.upperBound;
%                 g = zeros(1, nObj);
%                 g(n) = -sin(2*x)*(uBound - lBound);
% 
%             elseif obj.linked == 1
%                 g = cos(x)^2*thisJac(obj.crossLinkNo, :);
%                 uBound = p(obj.crossLinkNo);
%                 g(n) = -sin(2*x)*(uBound - lBound);
% 
%             elseif obj.linked == 2
%                 g = -cos(x)^2*thisJac(obj.crossLinkNo, :);
%                 uBound = obj.upperBound - p(obj.crossLinkNo);
%                 g(n) = -sin(2*x)*(uBound - lBound);
%             else
%                 error('Unrecognised value %d for linked property', obj.linked);
%             end
% 
% 
%         case 'squared'
%             g = zeros(1, nObj);
%             g(n) = 2*x;
% 
%         case 'linear'
%             g = zeros(1, nObj);
%             g(n) = 1;
% 
%         case 'dummy'
%             g = obj.gval(thisJac, p, obj.crossLinkNo);
% 
%         case 'constant'    
%             g = zeros(1, nObj);
%         otherwise
%             error('The linking type %s is not recognised.', obj.type);
%    end
% end % of grad