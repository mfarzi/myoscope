classdef cylinderGPD < compartment
    % CYLINDER 
    % (Cylinder Model with Gaussian Phase Distribution)
    %
    %   An CYLINDER object is a basic COMPARTMENT object represenitng 
    %   particles diffusion in a cylinder. This class encapsulates the 
    %   model parameters and basic methods for fitting the parameters to
    %   diffusion weigthed MR signals or synthesize signals for a given 
    %   diffusion scheme.
    %
    %   For mathematical background see
    %       Vangelderen, P., DesPres, D., Vanzijl, P.C.M. and Moonen, C.,
    %       "Evaluation of restricted diffusion in cylinders.
    %       Phosphocreatine in rabbit leg muscle.", Journal of Magnetic
    %       Resonance, Series B, 103(3), pp.255-260, 1994.
    %
    %   properties:
    %       name              - class name
    %       s0                - b0-signal with no diffusion weight or
    %                           signal attenuation at b-val=0 
    %       diffPar           - diffusivity along the cylinder axis [m^2/s]
    %       r                 - radius of ellipsoid along the major axis
    %                           [m]
    %       theta             - elevation angle; the angle between the
    %                           cylinder axis (n1) and the z-axis. [radian]
    %       phi               - azimuth angle; the angle between the x-axis
    %                           and the n1 projection onto the xy-plane. 
    %                           [radian]
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
    %       getCylinderAxis   - return the normal vector (n1) giving the
    %                           cylinder axis orientation.
    %
    %   See also: compartment, ellipticalCylinder, stick
    %
    % Mohsen Farzi
    % Email: m.farzi@leeds.ac.uk
    
    properties (SetAccess = 'protected')
        name   = 'cylinderGPD'; % class name
        s0     = 1;          % b0-signal with no diffusion weight

        diffPar   = 1e-9;       % diffusivity along the cylinder axis [m^2/s]

        r      = 5e-6;       % Radius of ellipsoid along the major axis [m]

        theta  = 0;          % elevation angle; the angle between the 
                             % cylinder axis (n1) and the z-axis. [radian]

        phi    = 0;          % azimuth angle; the angle between the x-axis
                             % and the n1 projection  onto the xy-plane. [radian]   
                             
        links     = [];     % vector of type LINKER that maps 
                            % constrained model parameters to
                            % unconstrained optimisation
                            % variables                     
    end
    
    properties (Access=protected)
        nParams = 5;           % number of model parameters
        modelParams  = [];          % model parameters 
                                    % [s0; diffPar; r; theta; phi]   
        hyperparams = [];                            
    end
    
    methods 
        function obj = cylinderGPD(varargin) 
            %CYLINDERGPD Construct Function.
            %   cylinderGPD() construct an object with default set of
            %   parameters
            %
            %   cylinderGPD(params)construct an object with given initial 
            %   params; [diffPar; r; theta; phi; alpha]
            %
            %   cylinderGPD(..., <Parameter>, <Value>, ...) allow
            %   passing model parameters as parameter-value pairs. If
            %   params is provided as first argument, the param-value pairs
            %   will be ignored.
            % 
            obj.modelParams = [obj.s0; obj.diffPar; obj.r; obj.theta; obj.phi];               
            %obj.fixedParams = false(obj.nParams, 1);
            
            obj.links = [linker('type', 'cos'   , 'lowerBound', 0, 'upperBound', 1, 'initRange', [0, 1]);  ...
                         linker('type', 'cos'   , 'lowerBound', 0.1e-9, 'upperBound',  3e-9, 'initRange', [1e-10, 3e-9]);  ...
                         linker('type', 'cos'   , 'lowerBound', 0.1e-6, 'upperBound', 20e-6, 'initRange', [1e-6, 20e-6]);  ...
                         linker('type', 'linear', 'initRange', [0, pi/2]);  ...
                         linker('type', 'linear', 'initRange', [0, 2*pi])];
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
                p.addParameter('diffPar'     , obj.diffPar       , @obj.validationFCN_diff);
                p.addParameter('r'        , obj.r          , @obj.validationFCN_r);
                p.addParameter('theta'    , obj.theta      , @obj.validationFCN_theta);
                p.addParameter('phi'      , obj.phi        , @obj.validationFCN_phi);  
                p.addParameter('fitter'   , obj.fitter     , @obj.validationFCN_fitter);
                
                p.parse(varargin{:});              
                
                % update the fitter
                obj.fitter = p.Results.fitter;
                     
                paramsIsProvided = not(ismember('params', p.UsingDefaults));
                
                if paramsIsProvided
                    obj.s0        = p.Results.params(1);
                    obj.diffPar      = p.Results.params(2);
                    obj.r         = p.Results.params(3);
                    obj.theta     = p.Results.params(4);
                    obj.phi       = p.Results.params(5);
                    
                    % warning if parameters are overwritten
                    paramList = {'s0', 'diffPar', 'r', 'theta', 'phi'};
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
                    obj.diffPar      = p.Results.diffPar;
                    obj.r         = p.Results.r;
                    obj.theta     = p.Results.theta;
                    obj.phi       = p.Results.phi;
                end
        
                obj.modelParams = [obj.s0; obj.diffPar; obj.r; ...
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
            
            obj.modelParams(2:end) = [obj.diffPar; obj.r; obj.theta; obj.phi];
        end % of rotateAxis
        
        
        function n = getCylinderAxis(obj)
            % compute the normal vector parallel to the cylinder axis
            n = [cos(obj.phi)*sin(obj.theta);...
                 sin(obj.phi)*sin(obj.theta);...
                 cos(obj.theta)];
        end
        %\\
    end % of methods (public)
        
    methods (Access = protected)        
        function s = synthesize(obj, scheme)
            % synthesize(obj, scheme) synthesize signals for diffusion
            % particles in acylinder. The signal is modeled as the product
            % of the signals parallel and perpendicular to the cylinder
            % axis
            %  
            G = [scheme.x, scheme.y, scheme.z].*scheme.G_mag;
            
            % compute the cylinderal axis
            n = obj.getCylinderAxis();
            
            % comput the gradient along and perpendicular
            % to the cylinder axis 
            Gpar_mag2 = (G*n).^2; 
            Gperp_mag2 = sum(G.^2,2)-Gpar_mag2;
            
            % compute the signal parallel to the cylinder axis
            Lpar = (scheme.DELTA-scheme.delta/3).*(obj.GAMMA*scheme.delta).^2*obj.diffPar;
            sPar = exp(-Lpar.*Gpar_mag2);

            % compute the signal perpendicular to the cylinder axis
            beta = obj.Jp1ROOTS/obj.r;          
            DB2 = obj.diffPar*beta.^2;
            nom = 2*scheme.delta*DB2' - 2  ...
                + 2*exp(-scheme.delta*DB2') ...
                + 2*exp(-scheme.DELTA*DB2') ...
                - exp(-(scheme.DELTA-scheme.delta)*DB2') ...
                - exp(-(scheme.DELTA+scheme.delta)*DB2');

            denom = obj.diffPar^2*beta.^6.*(obj.Jp1ROOTS.^2-1);
            
            Lperp = 2*obj.GAMMA^2*nom*(1./denom);
            
            sPerp = exp(-Lperp.*Gperp_mag2);
            
            % compute signal as the product of the value along and
            % perpendicular to the cylinder axis
            s = obj.s0* sPerp .* sPar;
            
        end % of synthesize
        
        function jac = jacobian(obj, scheme)
            % getJacobian return the gradient of singal wrt to optimisation
            % parameters, i.e. model parameters passed through the linkFun.
            %
            % NOTE: IF YOU CHANGE linkFun or invLinkFun, MAKE SURE TO EDIT
            % THIS CODE APPROPRIATELY AS WELL.
            %
            
            nScheme = size(scheme, 1);
            
            % initialise the jac with zeros
            jac = zeros(nScheme, obj.nParams);
            
%             x = obj.linkFun();                  % linking parameters           
            n = getCylinderAxis(obj);           % cylinder axis
            G = [scheme.x, scheme.y, scheme.z].*scheme.G_mag; % gradient vector
            Gpar_mag2 = (G*n).^2; 
            Gperp_mag2 = sum(G.^2,2)-Gpar_mag2;
            beta = obj.Jp1ROOTS/obj.r;          
            B2 = beta.^2;
            DB2 = obj.diffPar*beta.^2;
            
            Lpar = (scheme.DELTA-scheme.delta/3).*(obj.GAMMA*scheme.delta).^2*obj.diffPar;
           
            nom = 2*scheme.delta*DB2' - 2  ...
                + 2*exp(-scheme.delta*DB2') ...
                + 2*exp(-scheme.DELTA*DB2') ...
                - exp(-(scheme.DELTA-scheme.delta)*DB2') ...
                - exp(-(scheme.DELTA+scheme.delta)*DB2');

            denom = obj.diffPar^2*beta.^6.*(obj.Jp1ROOTS.^2-1);
            
            Lperp = 2*obj.GAMMA^2*nom*(1./denom);
            
            f0 = obj.s0*exp(-Lpar.*Gpar_mag2).*exp(-Lperp.*Gperp_mag2);
            
            % gradient wrt to s0
            jac(:,1) = f0/obj.s0;
            
            
            % gLpar_gDiff
            gLpar_gDiff = Lpar/obj.diffPar;
     
            % gLperp_gDiff
            gNom_gDiff = 2*scheme.delta*B2' ...
                       - 2*exp(-obj.diffPar*scheme.delta*B2').*(scheme.delta*B2') ...
                       - 2*exp(-obj.diffPar*scheme.DELTA*B2').*(scheme.DELTA*B2') ...
                       +   exp(-obj.diffPar*(scheme.DELTA-scheme.delta)*B2').*((scheme.DELTA-scheme.delta)*B2')...
                       +   exp(-obj.diffPar*(scheme.DELTA+scheme.delta)*B2').*((scheme.DELTA+scheme.delta)*B2');
            
            gDenom_gDiff = 2*obj.diffPar*beta.^6.*(obj.Jp1ROOTS.^2-1);
            tmp = bsxfun(@times, gNom_gDiff, denom') - bsxfun(@times, nom, gDenom_gDiff');
            tmp = bsxfun(@rdivide, tmp, denom'.^2);
            gLperp_gDiff = sum(tmp, 2)*2*obj.GAMMA^2;
            
            % gradient wrt diffPar      
            gf_gDiff = -(gLpar_gDiff.*Gpar_mag2 + gLperp_gDiff.*Gperp_mag2).*f0;
            jac(:,2) = gf_gDiff;
            
            % gLperp_gR
            gNom_gBm = 4*obj.diffPar*scheme.delta*beta' ...
                     - 2*exp(-obj.diffPar*scheme.delta*B2').*(2*obj.diffPar*scheme.delta*beta') ...
                     - 2*exp(-obj.diffPar*scheme.DELTA*B2').*(2*obj.diffPar*scheme.DELTA*beta') ...
                     +   exp(-obj.diffPar*(scheme.DELTA-scheme.delta)*B2').*(2*obj.diffPar*(scheme.DELTA-scheme.delta)*beta') ...
                     +   exp(-obj.diffPar*(scheme.DELTA+scheme.delta)*B2').*(2*obj.diffPar*(scheme.DELTA+scheme.delta)*beta');
            
            gDenom_gBm = 6*obj.diffPar^2*beta.^5.*(obj.Jp1ROOTS.^2-1);
            
            tmp = bsxfun(@times, gNom_gBm, denom') - bsxfun(@times, nom, gDenom_gBm');
            gLperp_gBm = 2*obj.GAMMA^2*bsxfun(@rdivide, tmp, denom'.^2);
            
            % gLperm_gR
            gBm_gR = -obj.Jp1ROOTS/obj.r^2;
            gLperp_gR = gLperp_gBm * gBm_gR;
            
            % gradient wrt R
            gf_gR = -(gLperp_gR.*Gperp_mag2).*f0;
            jac(:,3) = gf_gR;
            
            
            % Gradient wrt to theta
            gf_gN = (repmat(2*(Lperp-Lpar).*(G*n).*f0,1,3).*G);
            gN_gTheta = [ cos(obj.phi)*cos(obj.theta); ...
                          sin(obj.phi)*cos(obj.theta); ...
                         -sin(obj.theta)];
                     
            jac(:,4) = gf_gN*gN_gTheta;         
            
            % Gradient wrt to phi
            gN_gPhi   = [-sin(obj.phi)*sin(obj.theta); ...
                          cos(obj.phi)*sin(obj.theta); ...
                          0];
            jac(:,5) = gf_gN*gN_gPhi;                   
        end % of jacobian
        
        function updateParams(obj, p)
            obj.s0          = p(1);
            obj.diffPar        = p(2); % diffusivity along the cylinder axis [s/m^2]
            obj.r           = p(3); % Radius of the cylinder [m]
            obj.theta       = p(4); % the angle between cylinder axis n and the z-axis
            obj.phi         = p(5); % the angle between projection of n onto the xy-plan
                                    % and the x-axis 
            obj.modelParams = p;
        end    
    end % of methods (protected)
    
    methods (Static)
        function paramList = getParamsList()
            paramList = {'s0', 'diffPar', 'r', 'theta', 'phi'};
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
            assert(v(3) > 0, 'Value must be positive for "r".');
        end
        
        function validationFCN_s0(obj, v)
            assert(isnumeric(v) && isscalar(v) && (v > 0), ...
                   'Value must be scalar, numeric, and positive.');
        end
        
        function validationFCN_diff(obj, v)
            assert(isnumeric(v) && isscalar(v) && (v > 0), ...
                   'Value must be scalar, numeric, and positive.');
        end
        
        function validationFCN_r(obj, v)
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
end % of cylinder class

