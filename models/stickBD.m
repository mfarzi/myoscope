classdef stickBD < compartment
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
        name      = 'stickBD';     % class name
        s0        = 1;              % b0-signal with no diffusion weight
        
        diffPar   = 1e-9;           % diffusivity [m^2/s]
        
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
        
        kappa1 = 4;                 % Concentration paramter for Watson 
                                    % distribution
        
        kappa2 = 4;                 % Concentration paramter for Watson 
                                    % distribution                           
                                       
        links     = [];             % vector of type LINKER that maps 
                                    % constrained model parameters to
                                    % unconstrained optimisation
                                    % variables
    end
    
    properties (Access=private)
        % for faster computation of results
        cB = 0;                     % constant for bingmham distribution
        U  = zeros(3);              % unit frame
        gcB_gK1 = 0;
        gcB_gK2 = 0;
        ecDiffPar = 0;
        ecDiffPerp1 = 0;
        ecDiffPerp2 = 0;
    end
    
    properties (Access=protected)
        nParams = 7;           % number of model parameters
        modelParams  = [];     % model parameters 
                               % [s0; diffPar; diffPerp1; diffPerp2;
                               %  theta; phi; alpha; kappa1; kappa2]   
        hyperparams = [30; 30];                    
    end
    
    methods 
        function obj = stickBD(varargin) 
            %TENSOR Construct Function.
            %   stickBD() construct an object with default set of parameters
            %   DT = eye(3)*1e-9
            %
            %   stickBD(params) construct an object with given initial 
            %   params; [s0; diffPar; diffPerp1; diffPerp2; theta; phi; alpha]
            %
            %   stickBD(..., <Parameter>, <Value>, ...) allow
            %   passing model parameters as parameter-value pairs. If
            %   params is provided as first argument, the param-value pairs
            %   will be ignored.
            % 
            
            obj.updateParams([obj.s0; obj.diffPar; obj.theta; obj.phi; ...
                              obj.alpha; obj.kappa1; obj.kappa2]);               
            
            % define linkers for the zeppelin class
            obj.links = [linker('type', 'cos'   , 'lowerBound', 0, 'upperBound', 1, 'initRange', [0, 1]);  ...
                         linker('type', 'cos'   , 'lowerBound', 0.1e-9, 'upperBound', 3e-9, 'initRange', [1e-10, 3e-9]);  ...
                         linker('type', 'linear', 'initRange', [0, pi/2]);  ...
                         linker('type', 'linear', 'initRange', [0, 2*pi]);  ...
                         linker('type', 'linear', 'initRange', [0, pi]);    ...
                         linker('type', 'cos'    , 'lowerBound', 0, 'upperBound',  128, 'initRange', [1, 64]); ...
                         linker('type', 'cos'    , 'lowerBound', 0, 'upperBound',  128, 'initRange', [1, 64])];
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
                p.addParameter('theta'    , obj.theta      , @obj.validationFCN_theta);
                p.addParameter('phi'      , obj.phi        , @obj.validationFCN_phi);
                p.addParameter('alpha'    , obj.alpha      , @obj.validationFCN_alpha);   
                p.addParameter('fitter'   , obj.fitter     , @obj.validationFCN_fitter);
                p.addParameter('kappa1'   , obj.kappa1 );
                p.addParameter('kappa2'   , obj.kappa2 );
                p.addParameter('nBinsTheta', obj.hyperparams(1), @(v) v>0 && v == floor(v));
                p.addParameter('nBinsPhi'  , obj.hyperparams(2), @(v) v>0 && v == floor(v));
                
                p.parse(varargin{:});              
                
                % update the fitter
                obj.fitter = p.Results.fitter;
                obj.hyperparams(1) = p.Results.nBinsTheta;
                obj.hyperparams(2) = p.Results.nBinsPhi;
                
                paramsIsProvided = not(ismember('params', p.UsingDefaults));
                
                if paramsIsProvided
                    obj.s0        = p.Results.params(1);
                    obj.diffPar   = p.Results.params(2);
                    obj.theta     = p.Results.params(3);
                    obj.phi       = p.Results.params(4);
                    obj.alpha     = p.Results.params(5);
                    obj.kappa1    = p.Results.params(6);
                    obj.kappa2    = p.Results.params(7);
                    
                    % warning if parameters are overwritten
                    paramList = {'s0', 'diffPar', 'alpha'  , 'theta' , ...
                                 'phi', 'kappa1', 'kappa2'};
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
                    obj.theta     = p.Results.theta;
                    obj.phi       = p.Results.phi;
                    obj.alpha     = p.Results.alpha;
                    obj.kappa1    = p.Results.kappa1;
                    obj.kappa2    = p.Results.kappa2;
                end
        
                obj.updateParams([obj.s0; obj.diffPar; obj.theta; obj.phi; ... 
                                  obj.alpha; obj.kappa1; obj.kappa2]);             
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
            
            obj.theta = mod(obj.theta, 2*pi);
            obj.phi = mod(obj.phi, 2*pi);
            obj.alpha = mod(obj.alpha, 2*pi);
            
            if obj.theta > pi
                obj.theta = 2*pi - obj.theta;
                obj.phi = mod(pi+obj.phi, 2*pi);
                obj.alpha = mod(-obj.alpha, pi);
            end
            
            % make sure theta < pi/2
            if obj.theta>pi/2 
                obj.theta = pi - obj.theta;
                obj.phi = mod(pi + obj.phi, 2*pi);
                obj.alpha = mod(-obj.alpha, pi);
            end
            
            obj.modelParams(2:end) = [obj.diffPar; obj.theta; obj.phi; ...
                                      obj.alpha; obj.kappa1; obj.kappa2];
        end % of rotateAxis
         
        
        function U = getUnitFrame(obj)
            U = obj.U;
        end
        
        function DT = getDT(obj)
            DT = obj.U*diag([obj.ecDiffPar, obj.ecDiffPerp1, obj.ecDiffPerp2])*obj.U';
        end
        
        function setEcDiffParmas(obj)
            [probs, thisTheta, thisPhi] = discritiseBinghamDistribution(obj);    
            thisTheta = thisTheta(:);
            thisPhi = thisPhi(:);
            
            [x, y, z] = obj.sphericalToCanonical(thisTheta, thisPhi, 1);
            
            u1 = obj.U(:,1);
            u2 = obj.U(:,2);
            
            obj.ecDiffPar = sum((obj.diffPar*(u1(1)*x+u1(2)*y+u1(3)*z).^2).*probs.*sin(thisTheta), 1);
            obj.ecDiffPerp1 = sum((obj.diffPar*(u2(1)*x+u2(2)*y+u2(3)*z).^2).*probs.*sin(thisTheta), 1);
            obj.ecDiffPerp2 = obj.diffPar - obj.ecDiffPar - obj.ecDiffPerp1;
        end
        
        function [probs, thisTheta, thisPhi] = discritiseBinghamDistribution(obj)
            nBinsTheta = obj.hyperparams(1);
            nBinsPhi = obj.hyperparams(2);
            
            thetaVec = linspace(0, pi, nBinsTheta+1);
            phiVec = linspace(0, 2*pi, nBinsPhi+1);
            
            [thisTheta, thisPhi] = meshgrid((thetaVec(1:nBinsTheta)+thetaVec(2:nBinsTheta+1))/2,...
                (phiVec(1:nBinsPhi)+phiVec(2:nBinsPhi+1))/2);
            
            dS = (thetaVec(2:nBinsTheta+1)-thetaVec(1:nBinsTheta))'*...
                 (phiVec(2:nBinsPhi+1)-phiVec(1:nBinsPhi));
            dS = abs(dS/(4*pi));
            
            [x, y, z] = obj.sphericalToCanonical(thisTheta, thisPhi, 1);
            
            n = [x(:), y(:), z(:)];
            probs = obj.binghamPdf(n', obj.U(:,2), obj.U(:,3), -obj.kappa1, -obj.kappa2);
            probs = probs(:).*dS(:);
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
            
            s = obj.s0 * exp(-b.*(obj.ecDiffPar  *(G*obj.U(:,1)).^2 + ...
                                  obj.ecDiffPerp1*(G*obj.U(:,2)).^2 + ...
                                  obj.ecDiffPerp2*(G*obj.U(:,3)).^2));                
        end
        
        
        function jac = jacobian(obj, scheme)
            % jacobian(obj, scheme) return the gradient of signal wrt to
            % model parameters.
            %
            % see also: getJacobian
            %

            % initialise the jac with zeros
            jac = zeros(size(scheme,1), obj.nParams);
                       
            [probs, thisTheta, thisPhi] = discritiseBinghamDistribution(obj);    
            thisTheta = thisTheta(:);
            thisPhi = thisPhi(:);
            
            [x, y, z] = obj.sphericalToCanonical(thisTheta, thisPhi, 1);
            
            f0 = obj.synthesize(scheme);         % signal value 
            G = [scheme.x, scheme.y, scheme.z]; % gradient vector
            b = (scheme.DELTA-scheme.delta/3).* ...
                ((scheme.delta .*scheme.G_mag)*obj.GAMMA).^2; % b-value
                      
            % gradient wrt to s0
            jac(:,1) = f0/obj.s0;
            
            % gradient wrt dPar
            gf_gEcDiffPar = (-b.*((G*obj.U(:,1)).^2)).*f0;
            
            % gradient wrt dPerp1
            gf_gEcDiffPerp1 = (-b.*((G*obj.U(:,2)).^2)).*f0;
            
            % gradient wrt dPerp2
            gf_gEcDiffPerp2 = (-b.*((G*obj.U(:,3)).^2)).*f0;
            
            % gradient of dPar wrt diffPar
            gEcDiffPar_gDiffPar   = sum((obj.U(1,1)*x+obj.U(2,1)*y+obj.U(3,1)*z).^2.*probs.*sin(thisTheta), 1);
            gEcDiffPerp1_gDiffPar = sum((obj.U(1,2)*x+obj.U(2,2)*y+obj.U(3,2)*z).^2.*probs.*sin(thisTheta), 1);
            gEcDiffPerp2_gDiffPar = 1 - gEcDiffPar_gDiffPar - gEcDiffPerp1_gDiffPar;
            
            jac(:,2) = gf_gEcDiffPar.*gEcDiffPar_gDiffPar + gf_gEcDiffPerp1.*gEcDiffPerp1_gDiffPar + gf_gEcDiffPerp2.*gEcDiffPerp2_gDiffPar;
            
            % gradient wrt theta, phi, and alpha
            % using chain rule, comput gf_gV1, gf_gV2, gf_gV3.           
            gDPar_gU1 = sum(obj.diffPar*2*(obj.U(1,1)*x+obj.U(2,1)*y+obj.U(3,1)*z).*[x,y,z].*probs.*sin(thisTheta), 1);
            gDPar_gU2 = zeros(1,3);
            gDPar_gU3 = zeros(1,3);
            
            gDPerp1_gU1 = zeros(1,3);
            gP_gU2 = (-2*obj.kappa1)*probs.*(obj.U(1,2)*x+obj.U(2,2)*y+obj.U(3,2)*z).*[x,y,z];
            gP_gU3 = (-2*obj.kappa2)*probs.*(obj.U(1,3)*x+obj.U(2,3)*y+obj.U(3,3)*z).*[x,y,z];
            gDPerp1_gU2 = sum(obj.diffPar*2*(obj.U(1,2)*x+obj.U(2,2)*y+obj.U(3,2)*z).*[x,y,z].*probs.*sin(thisTheta), 1) + ...
                          sum((obj.diffPar*(obj.U(1,2)*x+obj.U(2,2)*y+obj.U(3,2)*z).^2).*sin(thisTheta).*gP_gU2, 1);
            gDPerp1_gU3 = sum((obj.diffPar*(obj.U(1,2)*x+obj.U(2,2)*y+obj.U(3,2)*z).^2).*sin(thisTheta).*gP_gU3, 1);
            
            gDPerp2_gU1 = - gDPar_gU1 - gDPerp1_gU1;
            gDPerp2_gU2 = - gDPar_gU2 - gDPerp1_gU2;
            gDPerp2_gU3 = - gDPar_gU3 - gDPerp1_gU3;
            
            gf_gU1 = repmat(-2*obj.ecDiffPar*b.*(G*obj.U(:,1)).*f0,1,3).*G   + ...
                     repmat(gf_gEcDiffPar, 1, 3).*gDPar_gU1          + ...
                     repmat(gf_gEcDiffPerp1, 1, 3).*gDPerp1_gU1      + ...
                     repmat(gf_gEcDiffPerp2, 1, 3).*gDPerp2_gU1;
                 
            gf_gU2 = repmat(-2*obj.ecDiffPerp1*b.*(G*obj.U(:,2)).*f0,1,3).*G + ...
                     repmat(gf_gEcDiffPar, 1, 3).*gDPar_gU2          + ...
                     repmat(gf_gEcDiffPerp1, 1, 3).*gDPerp1_gU2      + ...
                     repmat(gf_gEcDiffPerp2, 1, 3).*gDPerp2_gU2;
                 
            gf_gU3 = repmat(-2*obj.ecDiffPerp2*b.*(G*obj.U(:,3)).*f0,1,3).*G + ...
                     repmat(gf_gEcDiffPar, 1, 3).*gDPar_gU3          + ...
                     repmat(gf_gEcDiffPerp1, 1, 3).*gDPerp1_gU3      + ...
                     repmat(gf_gEcDiffPerp2, 1, 3).*gDPerp2_gU3;
            
            % gradient wrt theta
            gU1_gTheta = [ cos(obj.phi)*cos(obj.theta); ...
                           sin(obj.phi)*cos(obj.theta); ...
                          -sin(obj.theta)];
                      
            gU2_gTheta = [-cos(obj.alpha)*cos(obj.phi)*sin(obj.theta); ...
                          -cos(obj.alpha)*sin(obj.phi)*sin(obj.theta); ...
                          -cos(obj.alpha)*cos(obj.theta)];
                     
            gU3_gTheta = [ sin(obj.theta)*cos(obj.phi)*sin(obj.alpha);...
                           sin(obj.theta)*sin(obj.phi)*sin(obj.alpha);...
                           cos(obj.theta)*sin(obj.alpha)];
                       
            jac(:,3) = gf_gU1*gU1_gTheta + gf_gU2*gU2_gTheta + gf_gU3*gU3_gTheta;          
            
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
                     
            jac(:,4) = gf_gU1*gV1_gPhi   + gf_gU2*gV2_gPhi   + gf_gU3*gV3_gPhi;
            
            % gradient wrt alpha
            % gV1_gAlpha = zeros(3, 1);             
            
            gV2_gAlpha = [-sin(obj.alpha)*cos(obj.theta)*cos(obj.phi) - cos(obj.alpha)*sin(obj.phi);...
                          -sin(obj.alpha)*cos(obj.theta)*sin(obj.phi) + cos(obj.alpha)*cos(obj.phi);
                           sin(obj.alpha)*sin(obj.theta)];                     
           
            gV3_gAlpha = [-cos(obj.alpha)*cos(obj.theta)*cos(obj.phi)+sin(obj.alpha)*sin(obj.phi);...
                          -cos(obj.alpha)*cos(obj.theta)*sin(obj.phi)-sin(obj.alpha)*cos(obj.phi);...
                           cos(obj.alpha)*sin(obj.theta)];
           
            jac(:,5) = gf_gU2*gV2_gAlpha + gf_gU3*gV3_gAlpha;
            
            gP_gK1      = probs.*(obj.gcB_gK1/obj.cB-(obj.U(1,2)*x+obj.U(2,2)*y+obj.U(3,2)*z).^2);
            gP_gK2      = probs.*(obj.gcB_gK2/obj.cB-(obj.U(1,3)*x+obj.U(2,3)*y+obj.U(3,3)*z).^2);
            
            gDPar_gK1 = sum((obj.diffPar*(obj.U(1,1)*x+obj.U(2,1)*y+obj.U(3,1)*z).^2).*gP_gK1.*sin(thisTheta), 1);
            gDPerp1_gK1 = sum((obj.diffPar*(obj.U(1,2)*x+obj.U(2,2)*y+obj.U(3,2)*z).^2).*gP_gK1.*sin(thisTheta), 1);
            gDPerp2_gK1 = - gDPar_gK1 - gDPerp1_gK1;
            
            jac(:,6) = gf_gEcDiffPar.*gDPar_gK1 + gf_gEcDiffPerp1.*gDPerp1_gK1 + ...
                       gf_gEcDiffPerp2.*gDPerp2_gK1;
            
            gDPar_gK2 = sum((obj.diffPar*(obj.U(1,1)*x+obj.U(2,1)*y+obj.U(3,1)*z).^2).*gP_gK2.*sin(thisTheta), 1);
            gDPerp1_gK2 = sum((obj.diffPar*(obj.U(1,2)*x+obj.U(2,2)*y+obj.U(3,2)*z).^2).*gP_gK2.*sin(thisTheta), 1);
            gDPerp2_gK2 = - gDPar_gK2 - gDPerp1_gK2;
            
            jac(:,7) = gf_gEcDiffPar.*gDPar_gK2 + gf_gEcDiffPerp1.*gDPerp1_gK2 + ...
                       gf_gEcDiffPerp2.*gDPerp2_gK2;       
        end % of jacobian
        
        
            
        function updateParams(obj, p)
            obj.s0          = p(1);   % b0 signal
            obj.diffPar     = p(2);   % diffusivity [s/m^2]
            obj.theta       = p(3);   % elevation angle
            obj.phi         = p(4);   % azimuth angel
            obj.alpha       = p(5);   % angele beween v1_zrot and v2
            obj.kappa1      = p(6);
            obj.kappa2      = p(7);
            obj.modelParams = p;
            
            % update private variables
            obj.U = getUnitFrame(obj.theta, obj.phi, obj.alpha);
            obj.cB = obj.cBingham(-obj.kappa1, -obj.kappa2);
            [obj.gcB_gK1, obj.gcB_gK2] = obj.cBinghamGrad(-obj.kappa1, -obj.kappa2);
            obj.setEcDiffParmas();
            
        end    
        
        function updateHyperparams(obj, p)
            obj.hyperparams = p;
        end
        
        function p = binghamPdf(obj, n, u1, u2, kappa1, kappa2)
            p = exp(kappa1*(u1'*n).^2 + kappa2*(u2'*n).^2)/obj.cB;
        end
    end % of methods (protected)
    
    methods (Static)
        function paramList = getParamsList()
            paramList = {'s0', 'diffPar', 'theta', 'phi', 'alpha', ...
                         'kappa1', 'kappa2'};
        end
        
        function hyperparamsList = getHyperparamsList()
            hyperparamsList = {'nBinsTheta', 'nBinsPhi'};
        end
        
        function [x, y, z] = sphericalToCanonical(theta, phi, r)
            x = r.*sin(theta).*cos(phi);
            y = r.*sin(theta).*sin(phi);
            z = r.*cos(theta);
        end
        
        function cB = cBingham(k1, k2)
            cB = integral2(@(theta,phi) ...
                 exp(k1*sin(theta).^2.*sin(phi).^2+...
                     k2*cos(theta).^2).* ...
                     sin(theta),0,pi,0,2*pi)/(4*pi);
        end
        
        function [gcB_gK1, gcB_gK2] = cBinghamGrad(k1, k2)
            gcB_gK1 = integral2(@(theta,phi) ...
                                exp(k1*sin(theta).^2.*cos(phi).^2  + ...
                                k2*sin(theta).^2.*sin(phi).^2).* ...
                                sin(theta).^3.*cos(phi).^2, ...
                                0,pi,0,2*pi)/(4*pi);
                            
            gcB_gK2 = integral2(@(theta,phi) ...
                                exp(k1*sin(theta).^2.*cos(phi).^2  + ...
                                k2*sin(theta).^2.*sin(phi).^2).* ...
                                sin(theta).^3.*sin(phi).^2, ...
                                0,pi,0,2*pi)/(4*pi);                
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
