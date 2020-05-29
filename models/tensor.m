classdef tensor < compartment
    % TENSOR 
    %
    %   A TENSOR object is a basic COMPARTMENT object represenitng 
    %   particles diffusion in a medium using a diffusion tensor (DT). This
    %   class encapsulates six parameters to estimate the DT and basic
    %   methods for fitting the parameters to diffusion weigthed MR signals
    %   or synthesize signals for a given diffusion scheme.
    %
    %   For mathematical background see
    %       Basser, P.J., Mattiello, J. and LeBihan, D., "MR diffusion 
    %       tensor spectroscopy and imaging.", Biophysical journal, 66(1),
    %       pp.259-267, 1994.
    %
    %   properties:
    %       name              - class name
    %       s0                - b0-signal with no diffusion weight or
    %                           signal attenuation at b-val=0 
    %       diffPar           - parallel diffusivity along the primary axis
    %                           [m^2/s]
    %       diffPerp1         - perpendicular diffusivity along the
    %                           sceondary axis [m^2/s]
    %       diffPerp2         - perpendicualr diffusivity along the
    %                           tertiary axis [m^2/s]
    %       theta             - elevation angle; the angle between the
    %                           first eigen vector (v1) and the z-axis.
    %                           [radian]
    %       phi               - azimuth angle; the angle between the x-axis
    %                           and the v1 projection onto the xy-plane. 
    %                           [radian]
    %       alpha             - the angle between the second eigen vector
    %                           and the v1 rotated by pi/2 around the 
    %                           z-axis. [radian]
    %       fitter            - an "optimizer" object to fit model
    %                           parameters
    %
    %   methods (public):
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
    %       getADC            - return apparatent diffusion coefficient
    %       getFA             - return fractional anisotropy
    %       getEigenVec       - return three egien vectors [v1, v2, v3]
    %       getDT             - reurn the diffusion tensor   
    %
    %   See also: compartment, ball, zeppelin
    %
    % Mohsen Farzi
    % Email: m.farzi@leeds.ac.uk
    
    properties (SetAccess = 'protected')
        name      = 'tensor';       % class name
        s0        = 1;              % b0-signal with no diffusion weight
        
        diffPar   = 1e-9;           % diffusivity [m^2/s]
        
        diffPerp1 = 1e-9;           % diffusivity [m^2/s]
        
        diffPerp2 = 1e-9;           % diffusivity [m^2/s]
        
        theta     = 0;              % elevation angle; the angle between
                                    % the parimary axis (v1) and the
                                    % z-axis. [radian]                  
                                      
        phi       = 0;              % azimuth angle; the angle 
                                    % between the x-axis and the v1 
                                    % projection  onto the xy-plane.
                                    % [radian]

        alpha     = 0;              % the angle between the secondary 
                                    % eigen vector (v2) and the v1
                                    % rotation by pi/2 around the
                                    % xy-plane.     
                                    
        links     = [];             % vector of type LINKER that maps 
                                    % constrained model parameters to
                                    % unconstrained optimisation
                                    % variables
    end
    
    
    properties (Access=protected)
        nParams = 7;           % number of model parameters
        modelParams  = [];     % model parameters 
                               % [s0; diffPar; diffPerp1; diffPerp2;
                               %  theta; phi; alpha]   
        hyperparams = [];                       
    end
    
    methods 
        function obj = tensor(varargin) 
            %TENSOR Construct Function.
            %   tensor() construct an object with default set of parameters
            %   DT = eye(3)*1e-9
            %
            %   tensor(params) construct an object with given initial 
            %   params; [s0; diffPar; diffPerp1; diffPerp2; theta; phi; alpha]
            %
            %   tensor(..., <Parameter>, <Value>, ...) allow
            %   passing model parameters as parameter-value pairs. If
            %   params is provided as first argument, the param-value pairs
            %   will be ignored.
            % 
            
            obj.modelParams = [obj.s0; obj.diffPar; obj.diffPerp1; ...
                               obj.diffPerp2; obj.theta; obj.phi; ...
                               obj.alpha];               
            
            % define linkers for the zeppelin class
            obj.links = [linker('type', 'cos'   , 'lowerBound', 0, 'upperBound', 1, 'initRange', [0, 1]);  ...
                         linker('type', 'cos'   , 'lowerBound', 0.1e-9, 'upperBound', 3e-9, 'initRange', [1e-10, 3e-9]);  ...
                         linker('type', 'cos'   , 'lowerBound', 0.1e-9, 'upperBound', 3e-9, 'initRange', [1e-10, 3e-9]);  ...
                         linker('type', 'cos'   , 'lowerBound', 0.1e-9, 'upperBound', 3e-9, 'initRange', [1e-10, 3e-9]);  ...
                         linker('type', 'linear', 'initRange', [0, pi/2]);  ...
                         linker('type', 'linear', 'initRange', [0, 2*pi]);  ...
                         linker('type', 'linear', 'initRange', [0, pi])];
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
                p.addParameter('diffPerp1', obj.diffPerp1  , @obj.validationFCN_diffPerp1);
                p.addParameter('diffPerp2', obj.diffPerp2  , @obj.validationFCN_diffPerp2);
                p.addParameter('theta'    , obj.theta      , @obj.validationFCN_theta);
                p.addParameter('phi'      , obj.phi        , @obj.validationFCN_phi);
                p.addParameter('alpha'    , obj.alpha      , @obj.validationFCN_alpha);   
                p.addParameter('fitter'   , obj.fitter     , @obj.validationFCN_fitter);

                
                p.parse(varargin{:});              
                
                % update the fitter
                obj.fitter = p.Results.fitter;
                
                paramsIsProvided = not(ismember('params', p.UsingDefaults));
                
                if paramsIsProvided
                    obj.s0        = p.Results.params(1);
                    obj.diffPar   = p.Results.params(2);
                    obj.diffPerp1 = p.Results.params(3);
                    obj.diffPerp2 = p.Results.params(4);
                    obj.theta     = p.Results.params(5);
                    obj.phi       = p.Results.params(6);
                    obj.alpha     = p.Results.params(7);
                    
                    % warning if parameters are overwritten
                    paramList = {'s0'       , 'diffPar', 'diffPerp1', ...
                                 'diffPerp2', 'alpha'  , 'theta'    , ...
                                 'phi'};
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
                    obj.diffPerp1 = p.Results.diffPerp1;
                    obj.diffPerp2 = p.Results.diffPerp2;
                    obj.theta     = p.Results.theta;
                    obj.phi       = p.Results.phi;
                    obj.alpha     = p.Results.alpha;
                end
        
                obj.modelParams = [obj.s0; obj.diffPar; obj.diffPerp1; ...
                                   obj.diffPerp2; obj.theta; obj.phi; ... 
                                   obj.alpha];             
            end % of if nargin==0
        end % of set
        
        function rotateAxis(obj)
            % rotateAxis remove the ambiguity in estimation of parameters
            % theta, phi, and alpha by rotating the cooriante system of the
            % eigen vectors appropriately.
            % 
            % The orthonormal coordiante system V = [v1, v2, v3] should 
            % rotate such that the eigne vectors v1, v2, and v3 are 
            % associated with the eigen values in a descending order.
            %
            % More over, the direction of v1 or -v1 is selected such that
            % the angle theta is always in range [0, pi/2].
            % 
            % Finally, the direction v2 or -v2 is selected such that the
            % angle alpha is always between [0, phi].
            %
            % see also: getUnitFrame, getEulerAngles
            
            V = getUnitFrame(obj.theta, obj.phi, obj.alpha);
            
            % sort the eigen values in descending order
            eigVal = [obj.diffPar; obj.diffPerp1; obj.diffPerp2];
            [eigVal, sortID] = sort(eigVal, 'descend');
            
            obj.diffPar   = eigVal(1);
            obj.diffPerp1 = eigVal(2);
            obj.diffPerp2 = eigVal(3);
            
            V = V(:,sortID);
            % make sure V is still orthogonal following the permutation
            if uint8(norm(cross(V(:,1), V(:,2))-V(:,3))) ~= 0
                V(:,3) = -V(:,3);
            end
            [obj.theta, obj.phi, obj.alpha] = getEulerAngles(V);
            
            % make sure theta < pi/2
            if obj.theta>pi/2
                obj.theta = pi - obj.theta;
                obj.phi = pi + obj.phi;
                obj.alpha = -obj.alpha;
            end

            obj.phi = mod(obj.phi, 2*pi);
            
            % make sure alpha is in rage [0, pi]
            obj.alpha = mod(obj.alpha, pi);
            
            obj.modelParams(2:end) = [obj.diffPar; obj.diffPerp1;...
                                      obj.diffPerp2; obj.theta; ...
                                      obj.phi; obj.alpha];
        end % of rotateAxis
        
        
        function adc = getADC(obj)
            % compute Apparate Diffusion Coefficient
            adc = (obj.diffPar + obj.diffPerp1 + obj.diffPerp2)/3;
        end
        
        function fa = getFA(obj)
            % Fractional anisotropy
            lambda = [obj.diffPar, obj.diffPerp1, obj.diffPerp2];
            lambda_hat = mean(lambda);
            fa = sqrt(3/2)*norm(lambda-lambda_hat)/norm(lambda);
        end
                    
        function DT = getDT(obj)
            V = obj.getEigenVec();
            D = diag([obj.diffPar; obj.diffPerp1; obj.diffPerp2]);
            DT = V*D*V';
        end
        
        function V = getEigenVec(obj)
            V = getUnitFrame(obj.theta, obj.phi, obj.alpha);
        end
        %\\
    end % of method (public)
        
    methods (Access = protected)
        function s = synthesize(obj, scheme)
            % synthesize(obj, scheme) synthesize signals for diffusion
            % particles; s = s0 exp(-b*(G'*DT*G))
            % 
            G = [scheme.x, scheme.y, scheme.z];
            
            b = (scheme.DELTA-scheme.delta/3).* ...
                ((scheme.delta .*scheme.G_mag)*obj.GAMMA).^2;
            
            % compute the cylinderal axis
            V = obj.getEigenVec();
            
            s = obj.s0 * exp(-b.*(obj.diffPar  *(G*V(:,1)).^2 + ...
                                  obj.diffPerp1*(G*V(:,2)).^2 + ...
                                  obj.diffPerp2*(G*V(:,3)).^2)); 
        end
        
        function jac = jacobian(obj, scheme)
            % jacobian(obj, scheme) return the gradient of signal wrt to
            % model parameters.
            %
            % see also: getJacobian
            %

            % initialise the jac with zeros
            jac = zeros(size(scheme,1), obj.nParams);
            
            %x = obj.linkFun();                  % linking parameters           
            f0 = obj.synthesize(scheme);         % signal value 
            V = obj.getEigenVec();              % egine vectors
            G = [scheme.x, scheme.y, scheme.z]; % gradient vector
            b = (scheme.DELTA-scheme.delta/3).* ...
                ((scheme.delta .*scheme.G_mag)*obj.GAMMA).^2; % b-value
                      
            % gradient wrt to s0
            jac(:,1) = f0/obj.s0;
            
            % gradient wrt diffPar
            gf_gDiffPar = (-b.*((G*V(:,1)).^2)).*f0;
            jac(:,2) = gf_gDiffPar;
            
            % gradient wrt diffPerp1
            gf_gDiffPerp1 = (-b.*((G*V(:,2)).^2)).*f0;
            jac(:,3) = gf_gDiffPerp1;
            
            % gradient wrt diffPerp2
            gf_gDiffPerp2 = (-b.*((G*V(:,3)).^2)).*f0;
            jac(:,4) = gf_gDiffPerp2;
            
            % gradient wrt theta, phi, and alpha
            % using chain rule, comput gf_gV1, gf_gV2, gf_gV3.
            gf_gV1 = repmat(-2*obj.diffPar*b.*(G*V(:,1)).*f0,1,3).*G;
            gf_gV2 = repmat(-2*obj.diffPerp1*b.*(G*V(:,2)).*f0,1,3).*G;
            gf_gV3 = repmat(-2*obj.diffPerp2*b.*(G*V(:,3)).*f0,1,3).*G;
            
            % gradient wrt theta
            gV1_gTheta = [ cos(obj.phi)*cos(obj.theta); ...
                           sin(obj.phi)*cos(obj.theta); ...
                          -sin(obj.theta)];
                      
            gV2_gTheta = [-cos(obj.alpha)*cos(obj.phi)*sin(obj.theta); ...
                          -cos(obj.alpha)*sin(obj.phi)*sin(obj.theta); ...
                          -cos(obj.alpha)*cos(obj.theta)];
                     
            gV3_gTheta = [ sin(obj.theta)*cos(obj.phi)*sin(obj.alpha);...
                           sin(obj.theta)*sin(obj.phi)*sin(obj.alpha);...
                           cos(obj.theta)*sin(obj.alpha)];
                       
            jac(:,5) = gf_gV1*gV1_gTheta + gf_gV2*gV2_gTheta + gf_gV3*gV3_gTheta;          
            
            % gradient wrt phi
            gV1_gPhi   = [-sin(obj.phi)*sin(obj.theta); ...
                           cos(obj.phi)*sin(obj.theta); ...
                           0];
                       
            gV2_gPhi   = [-sin(obj.phi)*cos(obj.alpha)*cos(obj.theta) - cos(obj.phi)*sin(obj.alpha); ...          
                           cos(obj.phi)*cos(obj.alpha)*cos(obj.theta) - sin(obj.phi)*sin(obj.alpha); ...
                           0];
                       
            gV3_gPhi   = [ sin(obj.phi)*cos(obj.theta)*sin(obj.alpha) - cos(obj.phi)*cos(obj.alpha);...
                          -cos(obj.phi)*cos(obj.theta)*sin(obj.alpha) - sin(obj.phi)*cos(obj.alpha);...
                         0];    
                     
            jac(:,6) = gf_gV1*gV1_gPhi   + gf_gV2*gV2_gPhi   + gf_gV3*gV3_gPhi;
            
            % gradient wrt alpha
            % gV1_gAlpha = zeros(3, 1);             
            
            gV2_gAlpha = [-sin(obj.alpha)*cos(obj.theta)*cos(obj.phi) - cos(obj.alpha)*sin(obj.phi);...
                          -sin(obj.alpha)*cos(obj.theta)*sin(obj.phi) + cos(obj.alpha)*cos(obj.phi);
                           sin(obj.alpha)*sin(obj.theta)];                     
           
            gV3_gAlpha = [-cos(obj.alpha)*cos(obj.theta)*cos(obj.phi)+sin(obj.alpha)*sin(obj.phi);...
                          -cos(obj.alpha)*cos(obj.theta)*sin(obj.phi)-sin(obj.alpha)*cos(obj.phi);...
                           cos(obj.alpha)*sin(obj.theta)];
           
            jac(:,7) = gf_gV2*gV2_gAlpha + gf_gV3*gV3_gAlpha;
           
        end % of jacobian
        
            
        function updateParams(obj, p)
            obj.s0          = p(1);   % b0 signal
            obj.diffPar     = p(2);   % diffusivity [s/m^2]
            obj.diffPerp1   = p(3);   % diffusivity [s/m^2]
            obj.diffPerp2   = p(4);   % diffusivity [s/m^2]
            obj.theta       = p(5);   % elevation angle
            obj.phi         = p(6);   % azimuth angel
            obj.alpha       = p(7);   % angele beween v1_zrot and v2
            obj.modelParams = p;
        end    
        
        function updateHyperparams(obj, p)
            % do nothing
        end
    end % of methods (protected)
    
    methods (Static)
        function paramList = getParamsList()
            paramList = {'s0', 'diffPar', 'diffPerp1', ...
                         'diffPerp2', 'theta', 'phi', 'alpha'};
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
            assert(v(3) > 0, 'Value must be positive for "diffPerp1".');
            assert(v(4) > 0, 'Value must be positive for "diffPerp2".');
        end
        
        function validationFCN_s0(obj, v)
            assert(isnumeric(v) && isscalar(v) && (v > 0), ...
                   'Value must be scalar, numeric, and positive.');
        end
        
        function validationFCN_diffPar(obj, v)
            assert(isnumeric(v) && isscalar(v) && (v > 0), ...
                   'Value must be scalar, numeric, and positive.');
        end
        
        function validationFCN_diffPerp1(obj, v)
            assert(isnumeric(v) && isscalar(v) && (v > 0), ...
                   'Value must be scalar, numeric, and positive.');
        end
        
        function validationFCN_diffPerp2(obj, v)
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
        
        function validationFCN_alpha(obj, v)
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
end
