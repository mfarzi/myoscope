function p = getHyperparams(obj, hparamName)
% getHyperparams is a method for class COMPARTMENT
%
%   getHyperparams(obj) return the class property "hyperparams"
%
%   getHyperparams(obj, hparamName) return value for the specific 
%   hyperparameter. See getHyperparamsName for a full list of hyper-parameters.
%
% getHyperparams(obj)
if nargin == 1
    p = obj.hyperparams;
    return;
end

% getHyperparams(obj, hparamName)
assert(isa(hparamName, 'char'), 'MATLAB:compartment:getHyperparams',...
    'Input hyper-parameter name must be of type char');

% get model parameter names
paramsList = obj.getHyperparamsName;
idx = strcmp(paramsList, hparamName);
if any(idx)
    p = obj.hyperparams(idx);
else
    error('MATLAB:compartment:getHyperparams',...
        '"%s" is not a valid hyper-parameter name.', paramName);
end
end