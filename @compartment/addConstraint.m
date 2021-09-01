function addConstraint(obj, str)
    % addConstraint is a method for class MULTICOMPARTMENT
    %
    % addConstraint(obj, str) enforce the input constraint on model 
    % parameters. 
    %       input arguments:
    %
    %                   obj: Multicompartment object
    %                   str: Input of type char. Semantic parameter names must 
    %                        be used. See getParamsName(obj) for a list of
    %                        valid semantic names.
    %                        Examples:
    %                           1) 's0 = 1' set s0 as constant with value 1
    %                           2) 'diffPar >= 1.5e-9'
    %                           3) 'diffPar <= 1.5e-9'
    %                           4) 'diffPar >= diffPerp1'
    %                           5) 'diffPerp2 = diffPerp1'
    %
    % See also: getConstraints

    obj.links.addConstraint(str);
end