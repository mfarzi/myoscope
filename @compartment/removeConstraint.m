function removeConstraint(obj, varName)
    % removeConstraint is a method for class MULTICOMPARTMENT
    %
    %   removeConstraint(obj, varName) enforce the input variabe to be 
    %   independent from other variables.
    %
    % See also: addConstraint, getConstraints

    obj.links.removeConstraint(varName);
end