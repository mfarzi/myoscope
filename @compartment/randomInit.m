function p = randomInit(obj, seed)
% randomInit initialise the model parameters randomly using 
% plausible biophysical ranges
if nargin == 2
    % make sure random number generation is repeatable
    rng(seed);
end

x = obj.links.initx();
p = obj.links.map(x);
%obj.updateParams(p)
end