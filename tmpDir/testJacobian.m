function varargout = testJacobian(obj, scheme)
    jac_num = getJacobian(obj, scheme, true);

    jac_anl = getJacobian(obj, scheme, false);

    TOL = 1e-5;

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
        