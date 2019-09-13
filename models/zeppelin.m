classdef zeppelin < compartment
    % ZEPPELIN
    %
    %   A ZEPPELIN object is a basic COMPARTMENT object represenitng  
    %   particles diffusion in an anisotropic, but cylinderically symmetric
    %   medium. This class encapsulates the model parameters and basic
    %   methods for fitting the parameters to diffusion weigthed MR signals
    %   or synthesize signals for a given diffusion scheme.
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
    %       diffPar           - parallel diffusivity [m^2/s]
    %       diffPerp          - perpendicular diffusivity [m^2/s]
    %       theta             - elevation angle; the angle between the
    %                           first eigen vector (v1) and the z-axis.
    %                           [radian]
    %       phi               - azimuth angle; the angle between the x-axis
    %                           and the v1 projection onto the xy-plane. 
    %                           [radian]
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
    %   See also: compartment, tensor, ball
    %
    % Mohsen Farzi
    % Email: m.farzi@leeds.ac.uk
    
    properties (SetAccess = 'protected')
        name     = 'zeppelin';        % class name
        
        s0       = 1;                 % b0-signal with no diffusion 
                                      % weight
        
        diffPar  = 1e-9;              % parallel diffusivity [m^2/s]
        
        diffPerp = 0.5e-9;            % perpendicular diffusivity [m^2/s]
        
        theta    = 0;                 % elevation angle; the angle between 
                                      % the axis (n) and the z-axis.
                                      % [radian]
                                           
        phi      = 0;                 % azimuth angle; the angle between
                                      % the x-axis and the n projection 
                                      % onto the xy-plane. [radian]  
                                      
        links     = [];               % vector of type LINKER that maps 
                                      % constrained model parameters to
                                      % unconstrained optimisation
                                      % variables
    end
    
    properties (Access=protected)
        nParams = 5;           % number of model parameters
        modelParams = [];      % [s0; diffPar; diffPerp; theta; phi]
        hyperparams = [];
    end
    
    methods 
        function obj = zeppelin(varargin) 
            %ZEPPELIN Construct Function.
            %   zeppelin() construct an object with default set of parameters
            %   DT = (diffPar-diffPerp)*n*n' + diffPerp*eye(3)
            %
            %   zeppelin(params) construct an object with given initial 
            %   params; [s0; diffPar; diffPerp; theta; phi]
            %
            %   zeppelin(..., <Parameter>, <Value>, ...) allow
            %   passing model parameters as parameter-value pairs. If
            %   params is provided as first argument, the param-value pairs
            %   will be ignored.
            % 
            
            obj.modelParams = [obj.s0; obj.diffPar; obj.diffPerp; ...
                               obj.theta; obj.phi];               
            
            % define linkers for the zeppelin class
            obj.links = [linker('type', 'cos'   , 'lowerBound', 0, 'upperBound', 1, 'initRange', [0, 1]);  ...
                         linker('type', 'cos'   , 'lowerBound', 0.1e-9, 'upperBound', 3e-9, 'initRange', [1e-10, 3e-9]);  ...
                         linker('type', 'cos'   , 'lowerBound', 0.1e-9, 'upperBound', 3e-9, 'initRange', [1e-10, 3e-9]);  ...
                         linker('type', 'linear', 'initRange', [0, pi/2]);  ...
                         linker('type', 'linear', 'initRange', [0, 2*pi])];
            obj.links.addConstraint('p2>=p3');
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
            end
        end % of constructor    
        
        function set(obj, varargin)
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
                p.addParameter('diffPar'  , obj.diffPar    , @obj.validationFCN_diffPar);
                p.addParameter('diffPerp' , obj.diffPerp   , @obj.validationFCN_diffPerp);
                p.addParameter('theta'    , obj.theta      , @obj.validationFCN_theta);
                p.addParameter('phi'      , obj.phi        , @obj.validationFCN_phi);
                p.addParameter('fitter'   , obj.fitter     , @obj.validationFCN_fitter);
                
                p.parse(varargin{:});              
                
                % update the fitter
                obj.fitter = p.Results.fitter;
                
                paramsIsProvided = not(ismember('params', p.UsingDefaults));
                
                if paramsIsProvided
                    obj.s0        = p.Results.params(1);
                    obj.diffPar   = p.Results.params(2);
                    obj.diffPerp  = p.Results.params(3);
                    obj.theta     = p.Results.params(4);
                    obj.phi       = p.Results.params(5);
                    
                    % warning if parameters are overwritten
                    paramList = {'s0'       , 'diffPar', 'diffPerp', ...
                                 'theta'    , 'phi'};
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
                    obj.diffPar   = p.Results.diffPar;
                    obj.diffPerp  = p.Results.diffPerp;
                    obj.theta     = p.Results.theta;
                    obj.phi       = p.Results.phi;
                end
        
                obj.modelParams = [obj.s0; obj.diffPar; obj.diffPerp; ...
                                   obj.theta; obj.phi];             
            end % of if nargin==0
        end % of set
        
        
        function rotateAxis(obj)
            % rotateAxis remove the ambiguity in estimation of parameters
            % theta and phi.
            %
            % the direction of v1 or -v1 is selected such that the angle
            % theta is always in range [0, pi/2].
            %
            % see also: getUnitFrame, getEulerAngles
            
            obj.theta = mod(obj.theta, 2*pi);
            obj.phi = mod(obj.phi, 2*pi);
            
            if obj.theta > pi
                obj.theta = 2*pi - obj.theta;
                obj.phi = mod(pi+obj.phi, 2*pi);
            end
            
            % make sure theta < pi/2
            if obj.theta>pi/2 
                obj.theta = pi - obj.theta;
                obj.phi = mod(pi + obj.phi, 2*pi);
            end
            
            obj.modelParams(4:5) = [obj.theta; obj.phi];
        end % of rotateAxis
        %\\
        
    end % methods (public)
        
    methods (Access = protected)
        function s = synthesize(obj, scheme)
            % synthesize(obj, scheme) synthesize signals for diffusion
            % particles; s = s0 exp(-b*(G'*DT*G))
            %
            
            G = [scheme.x, scheme.y, scheme.z];
            
            b = (scheme.DELTA-scheme.delta/3).* ...
                ((scheme.delta .*scheme.G_mag)*obj.GAMMA).^2;
            
            % compute the cylinderal axis
            n = [cos(obj.phi)*sin(obj.theta);...
                 sin(obj.phi)*sin(obj.theta);...
                 cos(obj.theta)];
            
            s = obj.s0 * exp(-b.*((obj.diffPar-obj.diffPerp)*(G*n).^2+obj.diffPerp)); 
             
        end % of synthesize
        
        function jac = jacobian(obj, scheme)
            % getJacobian return the gradient of singal wrt to optimisation
            % parameters, i.e. model parameters passed through the linkFun.
            %
            % NOTE: IF YOU CHANGE linkFun or invLinkFun, MAKE SURE TO EDIT
            % THIS CODE APPROPRIATELY AS WELL.
            %
            
            % initialise the jac with zeros
            jac = zeros(size(scheme,1), obj.nParams);
                      
            f0 = obj.synthesize(scheme);         % signal value 
            V1 = [cos(obj.phi)*sin(obj.theta);...
                 sin(obj.phi)*sin(obj.theta);...
                 cos(obj.theta)];
            G = [scheme.x, scheme.y, scheme.z]; % gradient vector
            b = (scheme.DELTA-scheme.delta/3).* ...
                ((scheme.delta .*scheme.G_mag)*obj.GAMMA).^2; % b-value
                      
            % gradient wrt to s0
            jac(:,1) = f0/obj.s0;
            
            % gradient wrt diffPar
            gf_gDiffPar = (-b.*((G*V1).^2)).*f0;
            jac(:,2) = gf_gDiffPar;
            
            % gradient wrt diffPerp
            gf_gDiffPerp = (-b.*(-(G*V1).^2+1)).*f0;
            jac(:,3) = gf_gDiffPerp;
            
            
            % gradient wrt theta and phi
            % using chain rule, comput gf_gV1, gf_gV2, gf_gV3.
            gf_gV1 = repmat(-2*(obj.diffPar-obj.diffPerp)*b.*(G*V1).*f0,1,3).*G;
            
            % gradient wrt theta
            gV1_gTheta = [ cos(obj.phi)*cos(obj.theta); ...
                           sin(obj.phi)*cos(obj.theta); ...
                          -sin(obj.theta)];
                       
            jac(:,4) = gf_gV1*gV1_gTheta;          
            
            % gradient wrt phi
            gV1_gPhi   = [-sin(obj.phi)*sin(obj.theta); ...
                           cos(obj.phi)*sin(obj.theta); ...
                           0];
                     
            jac(:,5) = gf_gV1*gV1_gPhi;
           
        end % of jacobian
        
        function updateParams(obj, p)
            if nargin==1
                % do nothing
            else
                obj.s0          = p(1);
                obj.diffPar     = p(2);
                obj.diffPerp    = p(3);
                obj.theta       = p(4);
                obj.phi         = p(5);
                obj.modelParams = p;
            end
        end % of updateParams
        
    end % of methods (protected)
    
    methods (Static)
        function paramList = getParamsList()
            paramList = {'s0', 'diffPar', 'diffPerp', 'theta', 'phi'};
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
            assert(v(2) > 0, 'Value must be positive for "diffPar".');
            assert(v(3) > 0, 'Value must be positive for "diffPerp".');
        end
        
        function validationFCN_s0(obj, v)
            assert(isnumeric(v) && isscalar(v) && (v > 0), ...
                   'Value must be scalar, numeric, and positive.');
        end
        
        function validationFCN_diffPar(obj, v)
            assert(isnumeric(v) && isscalar(v) && (v > 0), ...
                   'Value must be scalar, numeric, and positive.');
        end
        
        function validationFCN_diffPerp(obj, v)
            assert(isnumeric(v) && isscalar(v) && (v > 0), ...
                   'Value must be scalar, numeric, and positive.');
        end
        
        function validationFCN_theta(obj, v)
            assert(isnumeric(v) && isscalar(v), ...
                   'Value must be scalar and numeric.');
        end
        
        function validationFCN_phi(obj, v)
            assert(isnumeric(v) && isscalar(v), ...
                   'Value must be scalar and numeric.');
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
    end % of methods (private)
    %\\
end % of ZEPPELIN class
