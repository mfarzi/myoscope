function jac = Fjac(obj, x, scheme, hparams)
% Fjac is a method for class MULTICOMPARTMENT
%
% Fjac(obj, x, varargin) returns the jacobian of singal with respect to 
% optimisation variables (x)

params = obj.links.map(x);
paramsJac = obj.jacobian(params, scheme, hparams);
linkingJac = obj.links.jacobian(x);
jac = paramsJac * linkingJac;
end % of fgrad