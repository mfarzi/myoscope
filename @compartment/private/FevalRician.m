function L = FevalRician(obj, x, data, scheme)
% Feval is a method for class MULTICOMPARTMENT
%
% Feval(obj, x, varargin) returns the difference between synthesized DW-MR
% signal  and the observed data. 
% f[m] = s_m - d_m

sigma = x(1)^2+0.0001;
xp = x(2:end);
params = obj.links.map(xp);
s = obj.synthesize(params, scheme);
L = ricianLogLikelihood(data, s, sigma);
end