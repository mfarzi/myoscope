classdef cylinderBD < compartment
    % CYLINDERGDR
    % (Cylinder Model with Watson Distribution of cylinder axis)
    %
    %   An CYLINDERGDR object is a basic COMPARTMENT object represenitng 
    %   particles diffusion in a media with various radii drawn from a 
    %   Gamma distribution. This class encapsulates the model parameters 
    %   and basic methods for fitting the parameters to diffusion weigthed
    %   MR signals or synthesize signals for a given diffusion scheme.
    %
    %   For mathematical background see
    %       Y. Assaf, R. Z. Freidlin, G. K. Rohde, and P. J. Basser,
    %       ?New modeling and experimental framework to characterize 
    %       hindered and restricted water diffusion in brain white 
    %       matter,? Magn. Reson. Med., vol. 52, no. 5, pp. 965?978, 2004.
    %
    %   properties:
    %       name              - class name
    %       s0                - b0-signal with no diffusion weight or
    %                           signal attenuation at b-val=0 
    %       diff              - diffusivity along the cylinder axis [m^2/s]
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
        name   = 'cylinderBD';          % class name
        s0     = 1;                     % b0-signal with no diffusion weight

        diffPar= 1e-9;                  % diffusivity along the cylinder
                                        % axis [m^2/s]
        
        r      = 10e-6;                 % cylinder radius                                    

        theta  = 0;                     % elevation angle; the angle between the 
                                        % cylinder axis (n1) and the z-axis. [radian]

        phi    = 0;                     % azimuth angle; the angle between the x-axis
                                        % and the n1 projection  onto the xy-plane. [radian]   
        
        alpha  = 0;
        
        kappa1 = 4;                    % Concentration paramter for Watson 
                                       % distribution
        
        kappa2 = 0;                    % Concentration paramter for Watson 
                                       % distribution                                
                                        
        links     = [];                % vector of type LINKER that maps 
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
        scheme = [];
    end
    
    properties (Access=protected)
        nParams = 8;                    % number of model parameters
        modelParams  = [];              % model parameters 
                                        % [s0; diffPar; r; theta; phi]                                    
        hyperparams = [50; 50];          % [nBinsTheta, nBinsPhi]
    end
    
    methods 
        function obj = cylinderBD(varargin) 
            %CYLINDERGDR Construct Function.
            %   cylinderBD() construct an object with default set of
            %   parameters
            %
            %   cylinderBD(params)construct an object with given initial 
            %   params; [diffPar; r; theta; phi; alpha]
            %
            %   cylinderBD(..., <Parameter>, <Value>, ...) allow
            %   passing model parameters as parameter-value pairs. If
            %   params is provided as first argument, the param-value pairs
            %   will be ignored.
            % 
            obj.modelParams = [obj.s0; obj.diffPar; obj.r; ...
                               obj.theta; obj.phi; obj.alpha; obj.kappa1; obj.kappa2];               
            
            obj.links = [linker('type', 'cos'   , 'lowerBound', 0, 'upperBound', 1, 'initRange', [0, 1]);  ...
                         linker('type', 'cos'    , 'lowerBound', 0.1e-9, 'upperBound',  3e-9, 'initRange', [1e-10, 3e-9]);  ...
                         linker('type', 'cos'    , 'lowerBound', 1e-6, 'upperBound', 20e-6, 'initRange', [5e-6, 15e-6]);  ...
                         linker('type', 'linear' , 'initRange', [0, pi/2]);  ...
                         linker('type', 'linear' , 'initRange', [0, 2*pi]); ...
                         linker('type', 'linear' , 'initRange', [0, pi]); ...
                         linker('type', 'cos'    , 'lowerBound', 0, 'upperBound',  128, 'initRange', [1, 64]); ...
                         linker('type', 'cos'    , 'lowerBound', 0, 'upperBound',  128, 'initRange', [1, 64])];
            obj.links.setName(obj.getParamsList);
            
            if nargin >0
                if ~isa(varargin{1}, 'char')
                    obj.set('params', varargin{:});
%                     if ismember('params', varargin{:})
%                         warning(['The first optional input is ignored', ...
%                                  ' since model parameters are set', ...
%                                  ' using value pair for "params".\n']);
%                     end
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
                p.addParameter('diffPar'     , obj.diffPar       , @obj.validationFCN_diffPar);
                p.addParameter('kappa1'    , obj.kappa1 );
                p.addParameter('kappa2'    , obj.kappa2 );
                p.addParameter('r'        , obj.r         , @obj.validationFCN_r);
                p.addParameter('theta'    , obj.theta      , @obj.validationFCN_theta);
                p.addParameter('phi'      , obj.phi        , @obj.validationFCN_phi); 
                p.addParameter('alpha'    , obj.alpha        , @obj.validationFCN_phi); 
                p.addParameter('fitter'   , obj.fitter     , @obj.validationFCN_fitter);
                p.addParameter('nBinsTheta', obj.hyperparams(1), @(v) v>0 && v == floor(v));
                p.addParameter('nBinsPhi' , obj.hyperparams(2), @(v) v>0 && v == floor(v));
                p.addParameter('scheme' , obj.scheme);
                
                p.parse(varargin{:});              
                
                % update the fitter
                obj.fitter = p.Results.fitter;
                
                obj.hyperparams(1) = p.Results.nBinsTheta;
                
                obj.hyperparams(2) = p.Results.nBinsPhi;
                
                obj.scheme = p.Results.scheme;
                
                paramsIsProvided = not(ismember('params', p.UsingDefaults));
                
                if paramsIsProvided
                    obj.s0        = p.Results.params(1);
                    obj.diffPar   = p.Results.params(2);
                    obj.r         = p.Results.params(3);
                    obj.theta     = p.Results.params(4);
                    obj.phi       = p.Results.params(5);
                    obj.alpha     = p.Results.params(6);
                    obj.kappa1    = p.Results.params(7);
                    obj.kappa2    = p.Results.params(8);
                    
                    % warning if parameters are overwritten
                    paramList = {'s0', 'diffPar', 'r', 'theta', 'phi', 'alpha', 'kappa1', 'kappa2'};
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
                    obj.r         = p.Results.r;
                    obj.theta     = p.Results.theta;
                    obj.phi       = p.Results.phi;
                    obj.alpha     = p.Results.alpha;
                    obj.kappa1    = p.Results.kappa1;
                    obj.kappa2    = p.Results.kappa2;
                end
        
                obj.updateParams([obj.s0; obj.diffPar; obj.r; ...
                                   obj.theta; obj.phi; obj.alpha; obj.kappa1; obj.kappa2]);             
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
            
            obj.modelParams(2:end) = [obj.diffPar; obj.r; obj.theta; obj.phi; obj.alpha; obj.kappa1; obj.kappa2];
        end % of rotateAxis
        
        
        function n = getCylinderAxis(obj)
            % compute the normal vector parallel to the cylinder axis
            n = [cos(obj.phi)*sin(obj.theta);...
                 sin(obj.phi)*sin(obj.theta);...
                 cos(obj.theta)];
        end
        
        function V = getUnitFrame(obj)
            V = getUnitFrame(obj.theta, obj.phi, obj.alpha);
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
    end % of methods (public)
        
    methods (Access = protected)        
        function [s, tmpAccessMemory] = synthesize(obj, scheme)
            % synthesize(obj, scheme) synthesize signals for diffusion
            % particles in acylinder. The signal is modeled as the product
            % of the signals parallel and perpendicular to the cylinder
            % axis
            %
            
            [probs, thisTheta, thisPhi] = discritiseBinghamDistribution(obj);     
            % compute cylinder axis for each (theta, phi) [3 x nInstance]
            n = [cos(thisPhi(:)).*sin(thisTheta(:)),...
                 sin(thisPhi(:)).*sin(thisTheta(:)),...
                 cos(thisTheta(:))]';
            
            % gradient matrix [nScheme x 3] for each measurment
            G_dir = [scheme.x, scheme.y, scheme.z];
            G_mag = scheme.G_mag;
            % comput the gradient along and perpendicular
            % to the cylinder axis [nScheme x nInstance]
            Gpar_mag = (G_dir*n).*G_mag;
            Gpar_mag2 = Gpar_mag.^2; 
            Gperp_mag2 = G_mag.^2-Gpar_mag2; 
            
            % compute the signal parallel to the cylinder axis
            % [nScheme x nInstance]
            Lpar = ((scheme.DELTA-scheme.delta/3).*(obj.GAMMA*scheme.delta).^2)*obj.diffPar;
            sPar = exp(-Lpar.*Gpar_mag2);
            
            % compute the signal perpendicular to the cylinder axis
            % [1 x length(Jp1ROOTS)]
            beta = obj.Jp1ROOTS'/obj.r;          
            B2   = beta.^2;
            
            
            %[nScheme x length(Jp1ROOTS)]
            delta_dot_beta2 = scheme.delta.*B2;
            DELTA_dot_beta2 = scheme.DELTA.*B2;
            
            nom = 2*obj.diffPar*delta_dot_beta2 - 2  ...
                + 2*exp(-obj.diffPar*delta_dot_beta2) ...
                + 2*exp(-obj.diffPar*DELTA_dot_beta2) ...
                - exp(-obj.diffPar*(DELTA_dot_beta2-delta_dot_beta2)) ...
                - exp(-obj.diffPar*(DELTA_dot_beta2+delta_dot_beta2));
            
            % [length(Jp1ROOTS) x 1]
            beta6_dot_r2beta2MinusOne = (beta.^6).*((obj.Jp1ROOTS').^2-1);
            denom = obj.diffPar^2 * beta6_dot_r2beta2MinusOne;
            
            %[nScheme x 1]
            Lperp = 2*obj.GAMMA^2*sum(nom./denom, 2);
            sPerp = exp(-Lperp.*Gperp_mag2);
            
            % compute signal as the product of the value along and
            % perpendicular to the cylinder axis
            sig = obj.s0 * sPerp .* sPar;
            
            % compute the average wrt Bingham Distribution
            s = sum(sig.*probs'.*sin(thisTheta(:)'), 2);
            
            if nargout == 2
                tmpAccessMemory.sig   = sig;
                tmpAccessMemory.Gpar_mag2 = Gpar_mag2;
                tmpAccessMemory.Gperp_mag2 = Gperp_mag2;
                
                tmpAccessMemory.delta_dot_beta2 = delta_dot_beta2;
                tmpAccessMemory.delta_dot_beta = scheme.delta.*beta;
                tmpAccessMemory.DELTA_dot_beta2 = DELTA_dot_beta2;
                tmpAccessMemory.DELTA_dot_beta = scheme.DELTA.*beta;
                tmpAccessMemory.beta6_dot_r2beta2MinusOne = beta6_dot_r2beta2MinusOne;
                tmpAccessMemory.beta5_dot_r2beta2MinusOne = (beta.^5).*((obj.Jp1ROOTS').^2-1);
                tmpAccessMemory.Lpar  = Lpar;
                tmpAccessMemory.Lperp = Lperp;
                tmpAccessMemory.nom   = nom;
                tmpAccessMemory.denom = denom;
                tmpAccessMemory.probs = probs;
                tmpAccessMemory.theta = thisTheta(:);
                tmpAccessMemory.phi   = thisPhi(:);
                tmpAccessMemory.n     = n;
            end
        end
        
        
        function jac = jacobian(obj, scheme)
            % getJacobian return the gradient of singal wrt to optimisation
            % parameters, i.e. model parameters passed through the linkFun.
            %
            % NOTE: IF YOU CHANGE linkFun or invLinkFun, MAKE SURE TO EDIT
            % THIS CODE APPROPRIATELY AS WELL.
            %
            [f0, tmpAccessMemory] = obj.synthesize(scheme);
            Lpar = tmpAccessMemory.Lpar;
            delta_dot_beta2 = tmpAccessMemory.delta_dot_beta2;
            delta_dot_beta  = tmpAccessMemory.delta_dot_beta;
            DELTA_dot_beta2 = tmpAccessMemory.DELTA_dot_beta2;
            DELTA_dot_beta = tmpAccessMemory.DELTA_dot_beta;
            beta6_dot_r2beta2MinusOne = tmpAccessMemory.beta6_dot_r2beta2MinusOne;
            beta5_dot_r2beta2MinusOne = tmpAccessMemory.beta5_dot_r2beta2MinusOne;
            
            nom = tmpAccessMemory.nom;
            denom = tmpAccessMemory.denom;
            
            Gpar_mag2 = tmpAccessMemory.Gpar_mag2;
            Gperp_mag2 = tmpAccessMemory.Gperp_mag2;
            
            sig = tmpAccessMemory.sig;
            
            probs = tmpAccessMemory.probs;
            thisTheta = tmpAccessMemory.theta;
            thisPhi = tmpAccessMemory.phi;
            %%
            % gSgS0 : gradient wrt s0
            jac(:,1) = f0/obj.s0;
            
            % gradient wrt diffPar
            % first compute gradients for each pair f(theta', phi') and then
            % average over the Bingham distriubtion.
            
            %gLpar_gDiffPar
            gLpar_gDiffPar = Lpar/obj.diffPar;
     
            % gLperp_gDiff
            gNom_gDiffPar = 2*delta_dot_beta2 ...
                       - 2*exp(-obj.diffPar*delta_dot_beta2).*delta_dot_beta2 ...
                       - 2*exp(-obj.diffPar*DELTA_dot_beta2).*DELTA_dot_beta2 ...
                       +   exp(-obj.diffPar*(DELTA_dot_beta2-delta_dot_beta2)).*(DELTA_dot_beta2-delta_dot_beta2)...
                       +   exp(-obj.diffPar*(DELTA_dot_beta2+delta_dot_beta2)).*(DELTA_dot_beta2+delta_dot_beta2);
            
            gDenom_gDiffPar = 2*obj.diffPar*beta6_dot_r2beta2MinusOne;
            tmp = (gNom_gDiffPar.*denom - nom.*gDenom_gDiffPar)./(denom.^2);
            gLperp_gDiff = 2*obj.GAMMA^2*sum(tmp, 2);
            
            % gradient wrt diffPar      
            gf_gDiffPar = -(gLpar_gDiffPar.*Gpar_mag2 + gLperp_gDiff.*Gperp_mag2).*sig;
            
            % gSgDiffPar: average gradient over the Bingham distribution
            jac(:,2) = sum(probs'.*sin(thisTheta').*gf_gDiffPar, 2);
            
            % gradient wrt r
            % first compute gradients for each pair f(theta', phi') and then
            % average over the Bingham distriubtion.
            % gLperp_gR
            gNom_gBm = 4*obj.diffPar*delta_dot_beta ...
                     - 2*exp(-obj.diffPar*delta_dot_beta2).*(2*obj.diffPar*delta_dot_beta) ...
                     - 2*exp(-obj.diffPar*DELTA_dot_beta2).*(2*obj.diffPar*DELTA_dot_beta) ...
                     +   exp(-obj.diffPar*(DELTA_dot_beta2-delta_dot_beta2)).*(2*obj.diffPar*(DELTA_dot_beta-delta_dot_beta)) ...
                     +   exp(-obj.diffPar*(DELTA_dot_beta2+delta_dot_beta2)).*(2*obj.diffPar*(DELTA_dot_beta+delta_dot_beta));
            
            gDenom_gBm = 6*obj.diffPar^2*beta5_dot_r2beta2MinusOne;
            
            gLperp_gBm = 2 * obj.GAMMA^2 * ...
                        (gNom_gBm.*denom - nom .* gDenom_gBm)./(denom.^2);
            
            % gLperp_gR
            gBm_gR = -(obj.Jp1ROOTS')/obj.r^2;
            gLperp_gR = sum(gLperp_gBm .* gBm_gR, 2);
            
            % gradient wrt R
            gf_gR = -(gLperp_gR.*Gperp_mag2).*sig;
            
            % gSgR
            jac(:,3)  = sum(probs'.*sin(thisTheta').*gf_gR, 2);
            
            %%          
            u2 = obj.U(:,2); u3 = obj.U(:,3);
            
            [nX, nY, nZ]  = obj.sphericalToCanonical(thisTheta, thisPhi, 1);
            u2DotN = u2(1)*nX + u2(2)*nY + u2(3)*nZ;
            u3DotN = u3(1)*nX + u3(2)*nY + u3(3)*nZ;
            
            gFb_gU2 = -2*obj.kappa1*u2DotN.*probs.*[nX, nY, nZ];
            gFb_gU3 = -2*obj.kappa2*u3DotN.*probs.*[nX, nY, nZ];
            
            % gradient wrt theta
            gU2_gTheta = [-cos(obj.alpha)*cos(obj.phi)*sin(obj.theta); ...
                          -cos(obj.alpha)*sin(obj.phi)*sin(obj.theta); ...
                          -cos(obj.alpha)*cos(obj.theta)];
                     
            gU3_gTheta = [ sin(obj.theta)*cos(obj.phi)*sin(obj.alpha);...
                           sin(obj.theta)*sin(obj.phi)*sin(obj.alpha);...
                           cos(obj.theta)*sin(obj.alpha)];
                       
            gFb_gTheta = gFb_gU2*gU2_gTheta + gFb_gU3*gU3_gTheta;   
            % gSgTheta 
            jac(:,4)   = sum(gFb_gTheta'.*sin(thisTheta').*sig, 2);
            
            % gradient wrt phi
            gU2_gPhi   = [-sin(obj.phi)*cos(obj.alpha)*cos(obj.theta) - cos(obj.phi)*sin(obj.alpha); ...          
                           cos(obj.phi)*cos(obj.alpha)*cos(obj.theta) - sin(obj.phi)*sin(obj.alpha); ...
                           0];
                       
            gU3_gPhi   = [ sin(obj.phi)*cos(obj.theta)*sin(obj.alpha) - cos(obj.phi)*cos(obj.alpha);...
                          -cos(obj.phi)*cos(obj.theta)*sin(obj.alpha) - sin(obj.phi)*cos(obj.alpha);...
                         0];    
                     
            gFb_gPhi = gFb_gU2*gU2_gPhi + gFb_gU3*gU3_gPhi;
            % gSgPhi
            jac(:,5)   = sum(gFb_gPhi'.*sin(thisTheta').*sig, 2);
            
            % gradient wrt alpha             
            gU2_gAlpha = [-sin(obj.alpha)*cos(obj.theta)*cos(obj.phi) - cos(obj.alpha)*sin(obj.phi);...
                          -sin(obj.alpha)*cos(obj.theta)*sin(obj.phi) + cos(obj.alpha)*cos(obj.phi);
                           sin(obj.alpha)*sin(obj.theta)];                     
           
            gU3_gAlpha = [-cos(obj.alpha)*cos(obj.theta)*cos(obj.phi)+sin(obj.alpha)*sin(obj.phi);...
                          -cos(obj.alpha)*cos(obj.theta)*sin(obj.phi)-sin(obj.alpha)*cos(obj.phi);...
                           cos(obj.alpha)*sin(obj.theta)];
                       
            gFb_gAlpha = gFb_gU2*gU2_gAlpha + gFb_gU3*gU3_gAlpha;           
            % gSgAlpha
            jac(:,6)   = sum(gFb_gAlpha'.*sin(thisTheta').*sig, 2);
            
            % gSgKappa1
            jac(:,7)   = sum(probs'.*sin(thisTheta').*sig.*(-u2DotN.^2+obj.gcB_gK1/obj.cB)', 2);
            
            % gSgKappa2
            jac(:,8)   = sum(probs'.*sin(thisTheta').*sig.*(-u3DotN.^2+obj.gcB_gK2/obj.cB)', 2);
        end
        
        
        function updateParams(obj, p)
            % update the primary model parameters
            obj.s0          = p(1);
            obj.diffPar     = p(2); % diffusivity along the cylinder axis [s/m^2]
            obj.r           = p(3); % Radius of the cylinder [m]
            obj.theta       = p(4); % the angle between cylinder axis n and the z-axis
            obj.phi         = p(5); % the angle between projection of n onto the xy-plan
                                    % and the x-axis 
            obj.alpha       = p(6);                        
            obj.kappa1      = p(7);  
            obj.kappa2      = p(8);  
            obj.modelParams = p;
            
            % update the intermediate model parameters (private variables)
            
            obj.U = getUnitFrame(obj.theta, obj.phi, obj.alpha);
            obj.cB = obj.cBingham(-obj.kappa1, -obj.kappa2);
            [obj.gcB_gK1, obj.gcB_gK2] = obj.cBinghamGrad(-obj.kappa1, -obj.kappa2);
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
            paramList = {'s0', 'diffPar', 'r', 'theta', 'phi', 'alpha', 'kappa1', 'kappa2'};
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
                 exp(k1*sin(theta).^2.*cos(phi).^2  + ...
                     k2*sin(theta).^2.*sin(phi).^2).* ...
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
        
        function g_kappa = hypergeomGrad(a, b, kappa)
            % M(a, b, kappa)= (gamma(b)/(gamma(b-a)*gamma(a))) ... 
            %           * integral(@(t) exp(kappa.*t).*t.^(a-1).*(1-t).^(b-a-1),0,1);
            g_kappa = (gamma(b)/(gamma(b-a)*gamma(a))) ...
                    * integral(@(t) exp(kappa.*t).*t.^(a).*(1-t).^(b-a-1),0,1);
        end
    end
    
    methods (Access = 'private')
        function validationFCN_params(obj, v)
            assert(numel(v) == obj.nParams && isnumeric(v),...
                   sprintf(['Value must be a numeric column or row',...
                            ' vector of size %d.\n'], obj.nParams));
                        
            assert(v(1) > 0, 'Value must be positive for "s0".');
            assert(v(2) > 0, 'Value must be positive for "diff".');
            assert(v(3) > 0, 'Value must be positive for "r".');
        end
        
        function validationFCN_s0(obj, v)
            assert(isnumeric(v) && isscalar(v) && (v > 0), ...
                   'Value must be scalar, numeric, and positive.');
        end
        
        function validationFCN_diffPar(obj, v)
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

