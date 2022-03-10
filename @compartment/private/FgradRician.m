function g = FgradRician(obj, x, data, scheme)
% Fjac is a method for class MULTICOMPARTMENT
%
% Fjac(obj, x, varargin) returns the jacobian of singal with respect to 
% optimisation variables (x)

%sigma = sigmaLink.map(x(1));
sigma = x(1)^2+0.0001;
xp = x(2:end);
params = obj.links.map(xp);

sig = obj.synthesize(params, scheme);
[gL_gSig, gL_gSigma] = ricianLogLikelihoodJacobian(data, sig, sigma);

gSig_gParams = obj.jacobian(params, scheme);
gParams_gX = obj.links.jacobian(xp);
gParams = (gSig_gParams*gParams_gX)'*gL_gSig;
g = [gL_gSigma*2*x(1); gParams];
end % of fgrad