function varargout = testJacobian(obj, params, schemefile)
    
    
    jac_anl = obj.jacobian(params, schemefile);
    
    % compute numerical jacobian
    jac_num = numericalJacobian(obj, params, schemefile);

    

    TOL = 1e-4;

    err = norm(jac_num(:)-jac_anl(:))/norm(jac_num(:));
    if err<TOL
        fprintf('Analytic solution is fine!\n');
    else
        fprintf('Analytic solution is buggy!\n');
    end

    if nargout == 0
        % do nothing!
    elseif nargout == 1
        varargout{1} = jac_anl;
    elseif nargout == 2
        varargout{1} = jac_anl;
        varargout{2} = jac_num;
    else
        error('Improper number of outputs.\n');
    end
end % of testjacobian

function jac = numericalJacobian(obj, params, schemefile)
    % getJacobian return the gradient of singal wrt to optimisation
    % parameters, i.e. model parameters passed through the linkFun.

    % compute the jacobian matrix using numerical technique.
    %\\
    EPS = 1e-4;                 % epsilon used for computation
                                % of numerical derivatives

    p0 = params;            % model params
    f0 = obj.synthesize(p0, schemefile);    % signal value
    
    % initialise the jac with zeros
    jac = zeros(schemefile.measurementsNum(), obj.nParams);
    
    for n = 1:obj.nParams
        p1 = p0;
        p1(n) = p0(n) + EPS*p0(n);
        jac(:,n) = (obj.synthesize(p1, schemefile) - f0)/(EPS*p0(n));
    end
    %\\
end %of getJacobian
        