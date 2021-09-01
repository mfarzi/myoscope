function g = grad(obj, X, P, jac)  
    % grad is a private method for the class linker
    % grad(obj, x, p, jac) compute the gradient of p[n] wrt x[1],...,x[N}
    %      Input Arguments:
    %                  obj: scalar LINKER object at index n
    %                    n: location of obj(n) in the vector
    %                    x: scalar unconstrained variable X[n]
    %                    P: constrained parameter vector
    %                  jac: current estimate of jacobian matrix [NxN]
    %
    %    Output Arguements:
    %                    g: gradient vector of size [1xN]
    %                       g = dP[n]/dX = [dP[n]/dx1,...dP[n]/dxN]
    
    N = obj.varNum;
    validateattributes(obj, {'linker'}, {'scalar'},...
        'MATLAB:linker:grad', 'obj');
    validateattributes(X, {'numeric'}, {'size', [N, 1]},...
        'MATLAB:linker:grad', 'X');
    validateattributes(P, {'numeric'}, {'size', [N, 1]},...
        'MATLAB:linker:grad', 'P');
    validateattributes(jac, {'numeric'}, {'size', [N, N]},...
        'MATLAB:linker:grad', 'jac');

    
     switch obj.type
        case 'independent'
            x = X(obj.index);
            u = obj.maximum;
            l = obj.minimum;
            
            g = zeros(1, N);
            g(obj.index) = obj.fgrad(x, u, l);
        case 'dependent'
            x = X(obj.index);
            u = obj.ubound(P(obj.dependence), obj.maximum, obj.minimum);
            l = obj.lbound(P(obj.dependence), obj.maximum, obj.minimum);
            
            gu= obj.uboundgrad(jac(obj.dependence,:), P(obj.dependence),obj.varNum);
            gl= obj.lboundgrad(jac(obj.dependence,:), P(obj.dependence),obj.varNum);
            
            g = zeros(1, N);
            g(obj.index) = obj.fgrad(x, u, l);
            g = g + cos(x)^2*(gu-gl) + gl;
        case 'dummy'
            g = obj.fgrad(jac(obj.dependence,:), P(obj.dependence));
        case 'constant'
            g = zeros(1, N);
        otherwise
            error('MATLAB:linker:scalarmap', ...
                strcat("Unknown type. Valid options:",...
                " 'free', 'dependent', 'dummy', or 'constant'."));
    end%of switch-case
end