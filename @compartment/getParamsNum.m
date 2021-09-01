function n = getParamsNum(obj, mode)
    % getParamsNum is a method for class MULTICOMPARTMENT
    %
    %   getParamsNum(obj) return the total number of model parameters
    %
    %   getParamsNum(obj, mode) 
    %           mode: String of type 'char': 
    %                                 'all': All parameters
    %                                'free': Optimisable variables
    %                            'constant': Not optimisable variables with
    %                                        constant values
    %                               'dummy': Not optimisable. An exact copy 
    %                                        of other variables.
    %                         'independent': Free variables which do not
    %                                        depend on other variables
    %                           'dependent': Free variables which depend on
    %                                        other variables
    %
    if nargin==1
        mode = 'all';
    end
    n = obj.links.getVarNum(mode);
end