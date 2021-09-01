classdef multicompartment < compartment 
    %MULTICOMPARTMENT 
    %
    %   A MULTICOMPARTMENT object allows combining one, two, or more basic 
    %   COMPARTMETN objects into a single model. This class encapsulates 
    %   the model parameters and basic methods for fitting the parameters
    %   to diffusion weigthed MR signals or synthesize signals for a given
    %   diffusion scheme.
    %
    %   For mathematical background see
    %       Panagiotaki, E., Schneider, T., Siow, B., Hall, M.G., Lythgoe,
    %       M.F. and Alexander, D.C., "Compartment models of the diffusion
    %       MR signal in brain white matter: a taxonomy and comparison.",
    %       Neuroimage, 59(3), pp.2241-2254, 2012.
    %
    %   properties (public):
    %       name           - Model name
    %       nParams        - Number of model parameters
    %       nHyperparams   - Number of hyper-parameters
    %
    %    properties (private):
    %       comp           - Cell array of constructing compartments
    %       links          - A column vector of class LINKER to map
    %                        unconstrained optimisation variables (x) to
    %                        constrained parameters (p).
    %                        p = obj.links.link(x);
    %                        x = obj.links.invLink(p);
    %
    %    methods:
    %       fit            - Fit model parameters to diffusion signals.
    %       synthesize     - Synthesize signals for the input diffusion scheme
    %       set            - Set model parameters or hyperparameters using
    %                        their names. Use getParamsList() to see the
    %                        semantic parameter names.
    %       randomInit     - randomly initialise model parameters.
    %       getCost        - return the Root Mean Square (RMS) error of
    %                        the fitted model. see class COMPARTMENT.
    %       getParams      - return the model parameters in a column
    %                        vector. see class COMPARTMENT.
    %       getParamsNum   - return the number of model parameters.
    %                        see class COMPARTMENT.
    %
    %   See also: compartment, tensor, ball, zeppelin, stick, cylinder,
    %   ellipticalCylinder
    %
    % Mohsen Farzi
    % Email: m.farzi@leeds.ac.uk
    
    properties (SetAccess = 'protected')
        name = '';               % Class name derived from its compartments
        nCompartments = [];      % number of basic COMPARTMETNT objects
        nParams = [];            % Number of parameters
        nHyperparams = [];       % Number of                                                  
    end
    
    properties (Access='protected') 
        comp = [];               % List of COMPARTMENT object 
        links     = [];          % vector of type LINKER that maps 
                                 % constrained model parameters to
                                 % unconstrained optimisation
                                 % variables 
        hyperparams = [];
        hyperparamsName = {};
    end
    
    methods
        function obj = multicompartment(varargin) 
            %MULTICOMPARTMETNMODEL Construct Function.
            %
            %   multicompartment(name) constructs a compartment model using
            %   the input model name. 
            %
            %   multicompartment(comp1, comp2,...) constructs a combined
            %   model with given compartments.
            % 
            
            % check for the number of input arguments
            assert(nargin>0, 'MATLAB:multicompartment:construct',...
                   'Too few input arguments.');
            
            % check consistency for multicompartment(name)
            modelName = varargin{1};
            if ischar(modelName)
                assert(nargin==1, 'MATLAB:multicompartment:construct',...
                       'Too many input arguments.');
                model = multicompartment.str2model(modelName);
                if isa(model, 'multicompartment')
                    obj = model;
                    return; 
                else
                    obj = multicompartment(model);
                    return;
                end
            end
            
            % check consistency for multicompartment(comp1, comp2,...)
            obj.nCompartments = nargin;
            obj.comp = varargin';
            % check if input compartments are valid inputs
            for i = 1:obj.nCompartments
                validateattributes(obj.comp{i}, ...
                    {'compartment', 'multicompartment'},...
                    {},'multicompartment:construction');
            end
            
            % set the object name
            nameArray = cellfun(@(c) c.name, obj.comp, 'UniformOutput', false);
            obj.name = strcat('[',strjoin(nameArray, '-'),']');
            
            % count the total number of paramters 
            nParamsArray = cellfun(@(c) c.nParams, obj.comp);
            obj.nParams = 1 + sum(nParamsArray);
            
            % count the total number of hyper-paramters
            nHyperparamsArray = cellfun(@(c) c.nHyperparams, obj.comp);
            obj.nHyperparams = sum(nHyperparamsArray);
            
            % vectorise all hyper-parameters into a single vector 
%             hyperparamsCell = cellfun(@(c) c.getHyperparams(),...
%                                  obj.comp, 'UniformOutput', false);
%             obj.hyperparams = cell2mat(hyperparamsCell);
            
            % update hyper-parameters name and values
            hparamsName = cell(obj.nHyperparams, 1);
            hparams = zeros(obj.nHyperparams, 1);
            iStart = 1;
            for i = 1:obj.nCompartments
                hparamsNum = obj.comp{i}.nHyperparams; 
                if hparamsNum>0
                    iEnd = iStart + hparamsNum -1;
                    [thisHparams,thisHparamsName] = obj.comp{i}.getHyperparams();
                    thisHparamsName = cellfun(@(x) sprintf('comp%d.%s', i, x),...
                        thisHparamsName, 'UniformOutput', false);
                    hparamsName(iStart:iEnd) = thisHparamsName;
                    hparams(iStart:iEnd) = thisHparams;
                    iStart = iEnd + 1;
                end
            end
            obj.hyperparamsName = hparamsName;
            obj.hyperparams = hparams;
            
            % vectorise all linking functions into a single vector
            linksArray = linker('s0', 'bounded', 0, 1);
            for i = 1:obj.nCompartments
                thisCompLinks = obj.comp{i}.getLinks();
                thisCompLinks(1).setName(strrep(thisCompLinks(1).getName(), 's0', 'vol'));
                thisCompParamsName = cellfun(@(c) sprintf('comp%d.%s', i, c), thisCompLinks.getName(), 'UniformOutput', false);
                thisCompLinks.setName(thisCompParamsName);
                linksArray = [linksArray; thisCompLinks];
            end
            obj.links = linksArray;
            
            % set constraint \sum_i comp_i = 1
            if obj.nCompartments==1
                obj.links.addConstraint('comp1.vol=1');
            elseif obj.nCompartments==2
                obj.links.addConstraint('comp2.vol=1-comp1.vol');
            elseif obj.nCompartments>2
                str = '1-comp1.vol';
                for n=2:obj.nCompartments-1
                    obj.links.addConstraint(sprintf('comp%d.vol<=%s', n, str));
                    str = strcat(str, sprintf('-comp%d.vol', n));
                end
                obj.links.addConstraint(sprintf('comp%d.vol=%s', obj.nCompartments, str));    
            end
        end % of constructor   
        
        function [sig, out] = synthesize(obj, params, scheme)
            % synthesize is a method for class MULTICOMPARTMENT
            %
            % synthesize(obj, scheme, params, hparams) synthesize DW-MR signal for the 
            % input diffusion scheme by adding the weighted signals from each 
            % compartment: s = s0 *sum_i f_i*s_i
            % 
            nScheme = size(scheme, 1);
            s = zeros(nScheme, 1);

            s0 = params(1);
            iStartParams = 2;
            for i = 1:obj.nCompartments
                iEndParams = iStartParams + obj.comp{i}.nParams()-1;

                s = s + obj.comp{i}.synthesize(params(iStartParams:iEndParams),...
                                               scheme);

                iStartParams  = iEndParams+1;
            end
            sig = s0*s;
            if nargout==2
                out.s = s;
            end
        end

        function jac = jacobian(obj, params, scheme)
            % jacobian is a method for class MULTICOMPARTMENT
            %
            % jacobian(obj,params, scheme, hparams) returns jacobian of synthesized
            % DW-MR signal with respect to model parameters. 
            %
            %       Input arguments:
            %                params: Numerical column vector of all model
            %                        parameters
            %                scheme: Diffusion scheme
            %
            %      output arguments:
            %                   jac: Numerical matrix of size M x nParams.
            %

            [f0, out] = obj.synthesize(params, scheme);  
            
            % read intermediate variables
            gf0_gs0 = out.s;

            % gradient wrt s0
            jac(:,1) = gf0_gs0;
            
            nScheme = size(scheme, 1);
            jac = zeros(nScheme, obj.nParams);
            
            s0  = params(1);
            iStartParams = 2;
            % gradient wrt to each compartment
            for i = 1: obj.nCompartments
                iEndParams = iStartParams + obj.comp{i}.nParams-1;

                gf_gThisComp = obj.comp{i}.jacobian(params(iStartParams:iEndParams),...
                                                    scheme);
                jac(:,iStartParams:iEndParams) = s0*gf_gThisComp;

                iStartParams = iEndParams + 1;
            end
        end%of jacobian
    end%of methods (public)
end%of class multicompartment
