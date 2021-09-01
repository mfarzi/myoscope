function p = getParams(obj, paramName)
% getParams is a method for class COMPARTMENT
%
%   getParams(obj) return the class property "params"
%
%   getParams(obj, paramName) return value for the specific parameter. See 
%   getParamsName for a full list of model parameters.
%

% getParams(obj)
if nargin == 1
    p = obj.params;
    return;
end

% getParams(obj, paramName)
assert(isa(paramName, 'char'), 'MATLAB:compartment:getParams',...
    'Input parameter name must be of type char');

% get model parameter names
paramsList = obj.getParamsName;
idx = strcmp(paramsList, paramName);
if any(idx)
    p = obj.params(idx);
else
    error('MATLAB:compartment:getParams',...
        '"%s" is not a valid parameter name.', paramName);
end
end%of getParams