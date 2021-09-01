function [hparams, hparamsList] = getHyperparams(obj)
% getParamsName is a method for class COMPARTMENT
%
%   getParamsName(obj) return the parameter names for the input object
%
hparams = obj.hyperparams;
hparamsList = obj.hyperparamsName();
end