function F = Feval(obj, x, data, scheme, hparams)
% Feval is a method for class MULTICOMPARTMENT
%
% Feval(obj, x, varargin) returns the difference between synthesized DW-MR
% signal  and the observed data. 
% f[m] = s_m - d_m

params = obj.links.map(x);
s = obj.synthesize(params, scheme, hparams);
F = s - data;
end