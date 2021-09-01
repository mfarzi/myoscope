function jac = getJacobian(obj, scheme, params, flag)
    % getJacobian return the gradient of singal wrt to optimisation
    % parameters, i.e. model parameters passed through the linkFun.

    if nargin == 3
        flag = false;
    end

    if flag
        % compute the jacobian matrix using numerical technique.
        %\\
        EPS = 1e-7;                 % epsilon used for computation
                                    % of numerical derivatives

        p0 = obj.getParams;            % model params
        x0 = obj.links.invLink(p0);    % linking parameters
        f0 = obj.synthesize(scheme);    % signal value

        nOptParams = length(x0);
        % initialise the jac with zeros
        jac = zeros(size(scheme,1), nOptParams);
        for n = 1:nOptParams
            x = x0; x(n) = x(n) + EPS;
            p = obj.links.link(x);
            obj.updateParams(p);
            jac(:,n) = (obj.synthesize(scheme) - f0)/EPS;
        end
        obj.updateParams(p0);
        %\\
    else
        % compute jacobian matrix using analytic solution
        paramsJac = obj.jacobian(scheme, params);
        x = obj.links.invLink(params);
        linkingJac = obj.links.jacobian(x);
        jac = paramsJac * linkingJac;
    end
    %\\
end %of getJacobian