function constraintList = getConstraints(obj)
% getConstraints is a method for class MULTICOMPARTMENT
%
%   getConstraints(obj) return a list of all constraints applied to model 
%   parameters. 
%
% See also: addConstraint

constraintList = obj.links.getConstraints();
end