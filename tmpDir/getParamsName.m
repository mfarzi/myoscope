function paramsList = getParamsName(obj)
% getParamsName is a method for class COMPARTMENT
%
%   getParamsName(obj) return the parameter names for the input object
%
paramsList = obj.links.getName();
end