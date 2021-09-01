function links = getLinks(obj)
% getLinks is a method for class COMPARTMENT
%
%   getLinks(obj) return a column vector of class LINKER used to establish
%   the relation between constrained model parameters and unconstrained 
%   optimisation parameters.
%
links = obj.links;
end