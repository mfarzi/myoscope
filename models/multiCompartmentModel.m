classdef multiCompartmentModel < compartment
    % TWOCOMPARTMETNMODEL (Two Compartment Model)
    %
    %   A TWOCOMPARTMETNMODEL object allows combining two basic 
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
    %   properties:
    %       name              - class name
    %       s0                - b0-signal with no diffusion weight or
    %                           signal attenuation at b-val=0 
    %       comp1             - First basic COMPARTMENT object
    %       comp2             - Second basic COMPARTMENT object
    %       vFractions        - Volume Fractions for each COMPARTMENT
    %       fitter            - an "optimizer" object to fit model
    %                           parameters
    %
    %    methods:
    %       fit               - fit model parameters to diffusion signals.
    %                           see the class COMPARTMENT.
    %       fitMultiRun       - fit model parameters using different set of
    %                           start points. see the class COMPARTMENT.
    %       synthesize        - synthesize signals for a diffusion scheme
    %                           see the class COMPARTMENT for definition.
    %       randomInit        - randomly initialise model parameters.
    %       set               - allow setting model parameters.
    %       fixParams         - Set specific parameters as constant during
    %                           model fitting. By default, all parameters
    %                           are set to be variables.
    %       getCost           - return the Root Mean Square (RMS) error of
    %                           the fitted model. see class COMPARTMENT.
    %       getParams         - return the model parameters in a column
    %                           vector. see class COMPARTMENT.
    %       getParamsNum      - return the number of model parameters.
    %                           see class COMPARTMENT.
    %       getFixedParams    - return a logical column vector stating if a
    %                           parameter is fixed (true) or not (false).
    %
    %   See also: compartment, tensor, ball, zeppelin, stick, cylinder,
    %   ellipticalCylinder
    %
    % Mohsen Farzi
    % Email: m.farzi@leeds.ac.uk
    
    properties (SetAccess = 'protected')
        name = '';               % class name derived from its compartments
        
        s0 = 1;                  % b0-signal with no diffusion weight
                                 % or volume fraction
        
        comp = [];               % List of COMPARTMENT object
        
        links     = [];          % vector of type LINKER that maps 
                                 % constrained model parameters to
                                 % unconstrained optimisation
                                 % variables
    end
    
    properties (Access=protected)
        nCompartments = [];      % number of basic COMPARTMETNT objects
        nParams = [];            % Number of parameters
        modelParams  = [];       % model parameters  
        hyperparams = [];
    end
    
    methods 
        function obj = multiCompartmentModel(varargin) 
            %TWOCOMPARTMETNMODEL Construct Function.
            %
            %   multiCompartmentModel(comp1, comp2) constructs a combined
            %   model with given compartments.
            %
            %   multiCompartmentModel(..., <Parameter>, <Value>, ...) 
            %   allows passing further information as parameter-value pairs.
            % 
            
            % set the number of compartments
            obj.nCompartments = nargin;
            
            % pre-allocate memory for each compartment
            compArray = cell(obj.nCompartments, 1);
            for i = 1:obj.nCompartments
                thisComp = varargin{i};
                if obj.isaValidCompartment(thisComp)
                    thisComp.set('s0', 1/nargin);
                    compArray{i,1} = thisComp;
                else
                    error('Not a valid compartment model.');
                end
            end
            obj.comp = compArray;
            
            % set the object name
            nameArray = cellfun(@(c) c.name, obj.comp, ...
                               'UniformOutput', false);
            obj.name = strcat('[',strjoin(nameArray, '-'),']');
            
            % count the total number of paramters
            paramsCellNum = cellfun(@(c) c.getParamsNum, obj.comp);
            obj.nParams = 1 + sum(paramsCellNum);
            
            % vectorise all parameters into a single vector 
            paramsCell = cellfun(@(c) c.getParams(), obj.comp, ...
                                 'UniformOutput', false);
            obj.modelParams = [obj.s0 ; cell2mat(paramsCell)];
            
            % vectorise all hyper-parameters into a single vector 
            paramsCell = cellfun(@(c) c.getHyperparams(),...
                                 obj.comp, 'UniformOutput', false);
            obj.hyperparams = cell2mat(paramsCell);
            
            % vectorise all linking functions into a single vector
            linkArray = cellfun(@(c) c.getLinks(), obj.comp, ...
                                'UniformOutput', false);
            tmpLink = linker('type',       'cos' , ...
                             'lowerBound',    0  , ...
                             'upperBound',    1  , ...
                             'initRange',  [0, 1]);           
            obj.links = [tmpLink; vertcat(linkArray{:})];
            
            obj.links.setName(obj.getParamsList);
            
            % set constraint \sum_i comp_i = 1
            str = sprintf('comp%d.s0+', 1:nargin);
            str(end) = '=';
            str = strcat(str, '1');
            obj.addConstraint(str);
        end % of constructor    
        
        function obj = set(obj, varargin)
            % set(obj, varargin) allow changing model parameters for an
            % existing object.
            %
            % obj.set() do nothing!
            %
            % obj.set(..., <Parameter>, <Value>, ...) 
            % allow passing model parameters as parameter-value pairs. 
            % Parse the parameter value pairs
            %
            if nargin==1
                % >> obj = obj.set()
            else
                p = inputParser;
                p.CaseSensitive = false;
            
                p.addParameter('s0'        , obj.s0         , @obj.validationFCN_s0); 
                p.addParameter('fitter'    , obj.fitter     , @obj.validationFCN_fitter);
                p.addParameter('params'    , obj.modelParams, @obj.validationFCN_modelParams); 
                
                p.parse(varargin{:});

                obj.fitter = p.Results.fitter;
                obj.s0 = p.Results.s0;
                
                obj.modelParams(1) = obj.s0;
                
                params = p.Results.params;
                if isrow(params)
                    params = params';
                end
                obj.updateParams(params);
            end % if (nargin == 1)
        end % of set
        
        function modelParamList = getParamsList(obj)
            modelParamList = {'s0'};
            for i = 1:obj.nCompartments
                thisCompParamList = obj.comp{i,1}.getParamsList();
                thisCompParamList = cellfun(@(x) sprintf('comp%d.%s', i, x), thisCompParamList, 'UniformOutput', false);
                
                modelParamList = [modelParamList, thisCompParamList];
            end
        end
        
        function hyperparamList = getHyperparamsList(obj)
            hyperparamList = {};
            for i = 1:obj.nCompartments
                thisCompHyperparamList = obj.comp{i,1}.getHyperparamsList();
                thisCompHyperparamList = cellfun(@(x) sprintf('comp%d.%s', i, x), thisCompHyperparamList, 'UniformOutput', false);
                
                hyperparamList = [hyperparamList, thisCompHyperparamList];
            end
        end
        
        function rotateAxis(obj)
            % rotateAxis remove the ambiguity in estimation of parameters
            % theta, phi, and alpha by rotating the cooriante system.
%             for i = 1:obj.nCompartments
%                 obj.comp{i}.rotateAxis;
%             end
            cellfun(@(c) c.rotateAxis(), obj.comp);
            
            % vectorise all parameters into a single vector 
            paramsCell = cellfun(@(c) c.getParams(), obj.comp, ...
                                 'UniformOutput', false);
            p = [obj.s0; cell2mat(paramsCell)];
            obj.updateParams(p);
        end
        
%         function p = getParams(obj, varargin)
%             IP = inputParser;
%             IP.CaseSensitive = false;
%             
%             validationFCN = @(v) assert(isscalar(v) && islogical(v), ...
%                 'Value must be scalar and logical.');
%             IP.addParameter('excludeS0', false , validationFCN);
% 
%             IP.parse(varargin{:});
%             
%             if IP.Results.excludeS0
%                 p = obj.modelParams;
%                 p(1) = [];
%             else
%                 p = obj.modelParams;
%             end 
%         end
        %\\
    end % of methods (public)
        
    methods (Access = protected)       
        function s = synthesize(obj, scheme)
            % synthesize(obj, scheme) synthesize signals for diffusion
            % particles in the compartment medium by adding the weighted 
            % signals from each compartment f_i*s_i
            % s = s0 *sum_i f_i*s_i
            % 
            nScheme = size(scheme, 1);
            s = zeros(nScheme, 1);
            for i = 1:obj.nCompartments
                s = s + obj.comp{i}.synthesize(scheme);
            end
            s = obj.s0*s;
        end
        
        function jac = jacobian(obj, scheme)
            %
            jac = zeros(size(scheme,1), obj.nParams);
            
            f0 = obj.synthesize(scheme);         % signal value
            
            % gradient wrt s0
            jac(:,1) = f0/obj.s0;
            
            idxStart = 1;
            % gradient wrt to each compartment
            for i = 1: obj.nCompartments
                gf_gThisComp = obj.s0*obj.comp{i}.getParamsJacobian(scheme);
                
                idxEnd = idxStart + obj.comp{i}.getParamsNum;
                idxStart = idxStart + 1;
                jac(:,idxStart:idxEnd) = gf_gThisComp;
                
                idxStart = idxEnd;
            end
            
        end % of jacobian
        
        function updateParams(obj, p)
            obj.s0 = p(1);
            
            idxStart = 1;
            for i = 1:obj.nCompartments
                idxEnd   = idxStart + obj.comp{i}.getParamsNum();
                idxStart = idxStart + 1;
             
                obj.comp{i}.updateParams(p(idxStart:idxEnd));
                idxStart = idxEnd;
            end
            
            obj.modelParams = p;
        end % of updateParams
        
        function updateHyperparams(obj, p)
            idxStart = 0;
            for i = 1:obj.nCompartments
                n = obj.comp{i}.getParamsNum('hyperparams');
                if n > 0
                    idxEnd   = idxStart + n;
                    idxStart = idxStart + 1;
             
                    obj.comp{i}.updateParams(p(idxStart:idxEnd));
                    idxStart = idxEnd;
                end
            end
            
            obj.hyperparams = p;
        end % of updateHyperparams
        %\\    
    end % of methods (protected)
    
    methods (Access = 'private')
        function state = isaValidCompartment(obj, v)
            state = isa(v, 'ellipticalCylinder') | ...
                    isa(v, 'cylinderGPD')           | ...    
                    isa(v, 'stick')                 | ...    
                    isa(v, 'ball')                  | ...
                    isa(v, 'tensor')                | ...
                    isa(v, 'zeppelin')              | ...
                    isa(v, 'cylinderGDR')           | ...
                    isa(v, 'cylinderWD')            | ...
                    isa(v, 'cylinderBD')            | ...
                    isa(v, 'tensorBD')              | ...
                    isa(v, 'stickBD')               | ...
                    isa(v, 'multiCompartmentModel');
        end
        
        function validationFCN_vFractions(obj, v)
            assert(numel(v) == obj.nCompartments && isnumeric(v) ...
                   && abs(1-sum(v)) < 1e4*eps && all(v>=0),...
                   sprintf(['Value must be a positive numeric column or row',...
                            ' vector of size %d with sum equals one.\n'], ...
                            obj.nCompartments));
        end
        
        function validationFCN_s0(obj, v)
            assert(isnumeric(v) && isscalar(v) && (v > 0), ...
                   'Value must be scalar, numeric, and positive.');
        end
        
        function validationFCN_fitter(obj, v)
            assert(isa(v, 'optimizer'), ...
                   'Value must be of type "optimizer".');
        end
        
        function validationFCN_fixParamsVec(obj, v)
            assert(numel(v) == obj.nCompartments && islogical(v),...
                   sprintf(['Value must be a logical column or row',...
                            ' vector of size %d.\n'], obj.nCompartments));
        end
        
        function validationFCN_modelParams(obj, v)
            isVector = iscolumn(v) || isrow(v);
            isCorrectSize = numel(v) == obj.nParams;
 
            if ~isVector || ~isCorrectSize || ~isnumeric(v)
                error('MATLAB:mulitCompartmentModel:set', ...
                     ['Value for "params" must be a double column or ',...
                      'row vector of size %d.\n'], obj.nParams);
            end
        end
        %\\
    end % of methods (private)
end
