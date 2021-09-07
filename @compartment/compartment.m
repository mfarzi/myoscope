classdef (Abstract) compartment < handle
    % COMPARTMENT 
    %
    %   COMPARTMENT is an ABSTRACT class providing a template definition 
    %   for methods and properties of a new single-compartment model. 
    %
    %   properties (TO BE IMPLEMENTED IN SUB-CLASSES):
    %       name               - Class name
    %       nParams            - Total number of model parameters
    %       nHyperparams       - Total number of hyper-parameters;
    %                            hyper-parameters are non-optimisable 
    %                            variables required for implementing
    %                            synthesize or jacobian methods.
    %
    %   methods (TO BE IMPLEMENTED IN SUB-CLASSES):
    %       synthesize         - Return DW-MR signal for the input scheme 
    %       jacobian           - Return the gradient of model wrt to its 
    %                            parameters. (analytic solution). 
    %       getLinks           - Return a column vector of class LINKER
    %                            that maps model parameters to unconstrained 
    %                            optimisation variables.
    %       getHyperparams     - Return hyper-parameters
    %   
    %   See also: tensor, ball, zeppelin, pancake, stick, cylinder, 
    %   ellipticalCylinder, multicompartment
    %
    % Mohsen Farzi
    % Email: m.farzi@leeds.ac.uk
    
    properties (Abstract, SetAccess = 'protected')  
        name;              % Given model name
        nCompartments;     % number of basic COMPARTMETNT objects
        nParams;           % Number of model parameters
        nHyperparams;      % Number of model hyper-parameters
    end
    
    properties (Abstract, Access = 'protected')
        comp;               % List of COMPARTMENT object 
        links;
        hyperparams;
        hyperparamsName;
    end
    
    methods (Abstract)
        s   = synthesize(obj, params, scheme, hyperparams);
        jac = jacobian(obj, params, scheme, hyperparams);
    end
    
    methods (Access=public)
        % utility
        links = getLinks(obj);
        paramsList = getParamsName(obj);
        [hparams, hparamsName] = getHyperparams(obj);
        setHyperparams(obj, varargin);
        n = getParamsNum(obj, mode);
        setParamsName(obj,varargin);
        
        writeModel(obj, filename, varargin);
        
        % optimisation methods       
        [params, rmse, exitFlag] = fit(obj, data, scheme, varargin);
        addConstraint(obj, str);
        constraintList = getConstraints(obj);
        p = randomInit(obj, seed);
        
        % test methods
        varargout = testJacobian(obj, params, scheme);
    end
    
    methods (Static)
        obj = readModel(filename);
        obj = str2model(modelName);
    end
end%of class