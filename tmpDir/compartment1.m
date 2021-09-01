classdef (Abstract) compartment
    % COMPARTMENT 
    %
    %   COMPARTMENT is an ABSTRACT class for defining new single or 
    %   multicompartment biophysical models. 
    %
    %   properties (TO BE IMPLEMENTED IN SUB-CLASSES):
    %       params             - Numerical column vector of all parameters
    %       nParams            - Total number of parameters
    %       hyperparams        - Numerical column vector of hyperparameters
    %       nHyperparams       - Total number of hyper-parameters
    %       hyperparamsName    - A list of hyper-parameters' names
    %       links              - A column vector of class LINKER to
    %                            establish the relation between constrained
    %                            model parameters and unconstrained 
    %                            optimisation parameters.
    %
    %   methods (TO BE IMPLEMENTED IN SUB-CLASSES):
    %       synthesize         - Return DW-MR signal for input scheme. 
    %       jacobian           - Return the gradient of model wrt to its 
    %                            parameters. (analytic solution). 
    %                            NOTE: A default numerical method is 
    %                                  alreadyimplemented in this class,
    %                                  but you are advised to overwrite
    %                                  this method.
    %
    %   methods (public)
    %       setParams          - Set model parameters
    %       getParams          - Get model parameters
    %       getParamsNum       - Return the number of model parameters
    %       getParamsName      - Return model parameters' names
    %       validateParams     - Assert input parameter complys with
    %                            constraints
    %       setHyperparams     - Set model hyper-parameters
    %       getHyperparams     - Get model hyper-parameters
    %       getHyperparamsNum  - Return the number of model hyper-parameters
    %       getHyperparamsName - Return hyper-parameters' names
    %       validateHyperparams- Assert hyper-parameters comply with
    %                            constraitns
    %       getLinks           - Return model links
    %
    %   
    %   See also: tensor, ball, zeppelin, cylinder, ellipticalCylinder,
    %   stick, multicompartment
    %
    % Mohsen Farzi
    % Email: m.farzi@leeds.ac.uk
    
    properties (Abstract, Constant)  
        name;              % Given model name
        %params;            % 
        nParams;           % Number of model parameters
        %hyperparams;       % Numerical column vector of hyperparameters  
        %hyperparamsName;   % List of hyper-param Names
        nHyperparams;      % Number of model hyper-parameters
    end
    
    methods (Abstract, Static)
        s   = synthesize(scheme, params, hyperparams);
        jac = jacobian(scheme, params, hyperparams);
        links = getLinks();
        paramsList = getParamsName();
        hparamsList = getHyperparamsName()
    end
    
    methods (Static)
        validateParams(obj, p);
        validateHyperparams(obj, p);
        
    end % of methods (public)
end%of class