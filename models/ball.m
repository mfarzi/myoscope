classdef ball < compartment
    % BALL
    %
    %   A BALL object is a basic COMPARTMENT object represenitng particles 
    %   diffusion in an isotropic, unconstrained medium. This class
    %   encapsulates the model parameters and basic methods for fitting  
    %   the parameters to diffusion weigthed MR signals or synthesize 
    %   signals for a given diffusion scheme.
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
    %       diff              - diffusivity of ball [m^2/s]
    %       fitter            - an "optimizer" object to fit model
    %                           parameters
    %
    %   methods:
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
    %   See also: compartment, tensor, zeppelin
    %
    % Mohsen Farzi
    % Email: m.farzi@leeds.ac.uk
    
    properties (SetAccess = 'protected')
        name = 'ball';              % class name
        
        s0 = 1;                     % b0-signal with no diffusion weight
        
        diff = 2.3e-9;              % diffusivity [m^2/s] 
        
        links     = [];             % vector of type LINKER that maps 
                                    % constrained model parameters to
                                    % unconstrained optimisation
                                    % variables
    end
    
    
    properties (Access=protected)
        nParams = 2;           % number of model parameters
        modelParams = [];           % [s0; diff]
        hyperparams = [];
    end
    
    methods 
        function obj = ball(varargin) 
            %BALL Construct Function.
            %   ball() construct an object with default set of parameters
            %   DT = eye(3)*2.3e-9
            %
            %   tensor(params) construct an object with given initial 
            %   params; [s0; diff]
            %
            %   tensor(..., <Parameter>, <Value>, ...) allow
            %   passing model parameters as parameter-value pairs. If
            %   params is provided as first argument, the param-value pairs
            %   will be ignored.
            % 
            obj.modelParams = [obj.s0; obj.diff];               
            
            % define linkers for the zeppelin class
            obj.links = [linker('type', 'cos'   , 'lowerBound', 0, 'upperBound', 1, 'initRange', [0, 1]);  ...
                         linker('type', 'cos'   , 'lowerBound', 0.1e-9, 'upperBound', 3e-9, 'initRange', [1e-10, 3e-9])];
            obj.links.setName(obj.getParamsList);
            
            if nargin >0
                if ~isa(varargin{1}, 'char')
                    obj.set('params', varargin{:});
                    if ismember('params', varargin{:})
                        warning(['The first optional input is ignored', ...
                                 ' since model parameters are set', ...
                                 ' using value pair for "params".\n']);
                    end
                else
                    obj.set(varargin{:});
                end
            end % of (nargin>0)
        end % of constructor    
        
        function obj = set(obj, varargin)
            % set(obj, varargin) allow changing model parameters for an
            % existing object.
            %
            % obj.set() do nothing!
            %
            % obj.set(..., <Parameter>, <Value>, ...) 
            % allow passing model parameters as parameter-value pairs. 
            %
            % NOTE: If 'params' are provided, the param-value pairs will be
            % ignored.
            %
            if nargin==1
                % >> obj = obj.set()
                % Do Nothing!
            else
                % Parse the parameter-value pairs
                p = inputParser;
                p.CaseSensitive = false;
                p.addParameter('params'   , obj.modelParams, @obj.validationFCN_params);
                p.addParameter('s0'       , obj.s0         , @obj.validationFCN_s0);                        
                p.addParameter('diff'     , obj.diff       , @obj.validationFCN_diff);
                p.addParameter('fitter'   , obj.fitter     , @obj.validationFCN_fitter);
                
                p.parse(varargin{:});              
                
                % update the fitter
                obj.fitter = p.Results.fitter;
                
                paramsIsProvided = not(ismember('params', p.UsingDefaults));
                
                if paramsIsProvided
                    obj.s0        = p.Results.params(1);
                    obj.diff   = p.Results.params(2);
                    
                    % warning if parameters are overwritten
                    paramList = {'s0', 'diff'};
                    id = find(~(ismember(paramList, p.UsingDefaults)));
                    if length(id) == 1
                        warning(['Parameters are already provided using ', ...
                                 'vector format. The input for "%s" ', ...
                                 'is ignored.'], paramList{id});
                    elseif length(id) > 1
                        str = sprintf('"%s"', paramList{id(1)});
                        for i = 2:length(id)
                            if i == length(id)
                                str = strcat(str, sprintf(' and "%s"', paramList{id(i)}));
                            else
                                str = strcat(str, sprintf(', "%s"', paramList{id(i)}));
                            end
                        end
                        warning(['Parameters are already provided using ', ...
                                 'vector format. The inputs for %s ', ...
                                 'are ignored.'], str);
                    end
                %\\    
                else % use param-value pairs to set the parameters
                    obj.s0        = p.Results.s0;
                    obj.diff      = p.Results.diff;
                end
        
                obj.modelParams = [obj.s0; obj.diff];
            end % of if nargin==0
            %\\
        end % of set(obj, varargin)
        
        function rotateAxis(obj)
            % rotateAxis remove the ambiguity in estimation of parameters
            % theta, phi, and alpha by rotating the cooriante system.
            %
            % This function is only for consitentcy accross compartment 
            % models!
            %
            % >> obj.roateAxis()
            % >> Do Nothing! 
        end
        %\\
    end % methods (public)
        
    methods (Access = protected)
        function s = synthesize(obj, scheme)
            % synthesize(obj, scheme) synthesize signals for diffusion
            % particles in an isotropic, unconstrained environment.
            %
            b = (scheme.DELTA-scheme.delta/3).* ...
                ((scheme.delta .*scheme.G_mag)*obj.GAMMA).^2;
            s = obj.s0*exp(-b*obj.diff);
        end % of synthesize
        
        function jac = jacobian(obj, scheme)
            % getJacobian return the gradient of singal wrt to optimisation
            % parameters, i.e. model parameters passed through the linkFun.
            %
            % NOTE: IF YOU CHANGE linkFun or invLinkFun, MAKE SURE TO EDIT
            % THIS CODE APPROPRIATELY AS WELL.
            %
            
            nScheme = size(scheme, 1);
            
            b = (scheme.DELTA-scheme.delta/3).* ...
                ((scheme.delta .*scheme.G_mag)*obj.GAMMA).^2;
            
            f0 = obj.s0*exp(-b*obj.diff);
            % initialise the jac with zeros
            jac = zeros(nScheme, obj.nParams);
            
            % gradient wrt s0
            jac(:,1) = exp(-b*obj.diff);
            
            % gradient wrt diff  
            gf_gDiff = -b.*f0;
            jac(:,2) = gf_gDiff;
            %\\
        end % of getParamsJacobian
        
        function updateParams(obj, p)
                obj.s0          = p(1);
                obj.diff        = p(2);
                obj.modelParams = p;
        end % of updateParams
        
    end % of methods (protected)
    
    methods (Static)
        function paramList = getParamsList()
            paramList = {'s0', 'diff'};
        end
        
        function hyperparamsList = getHyperparamsList()
            hyperparamsList = {};
        end
    end
    
     methods (Access = 'private')
        function validationFCN_params(obj, v)
            assert(numel(v) == obj.nParams && isnumeric(v),...
                   sprintf(['Value must be a numeric column or row',...
                            ' vector of size %d.\n'], obj.nParams));
                        
            assert(v(1) > 0, 'Value must be positive for "s0".');
            assert(v(2) > 0, 'Value must be positive for "diff".');
        end
        
        function validationFCN_s0(obj, v)
            assert(isnumeric(v) && isscalar(v) && (v > 0), ...
                   'Value must be scalar, numeric, and positive.');
        end
        
        function validationFCN_diff(obj, v)
            assert(isnumeric(v) && isscalar(v) && (v > 0), ...
                   'Value must be scalar, numeric, and positive.');
        end
            
        function validationFCN_fitter(obj, v)
            assert(isa(v, 'optimizer'), ...
                   'Value must be of type "optimizer".');
        end
        
        function validationFCN_fixParamsVec(obj, v)
            assert(numel(v) == obj.nParams && islogical(v),...
                   sprintf(['Value must be a logical column or row',...
                            ' vector of size %d.\n'], obj.nParams));
        end
        
        function validationFCN_fixParams(obj, v)
            assert(isscalar(v) && islogical(v),...
                   'Value must be scalar and logical.');
        end
     end % of method (private)
end % of BALL class
