classdef zeppelinBDA < compartment
    %ZEPPELINBDA 
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
    
    properties (SetAccess='protected')
        name = 'tensorBDA';    
        nCompartments = 1;     % number of basic COMPARTMETNT objects 
        nParams = 8;
        nHyperparams = 2;
    end
    
    properties (Access = 'protected')
        comp = [];             % List of COMPARTMENT object 
        links;
        hyperparams = [];
        hyperparamsName = {};
    end
    
    methods 
        function obj = zeppelinBDA() 
            %ZEPPELINBDA Construct Function.
            %
            %   zeppelinBDA() constructs a single-compartment model
            
            % initialise the linking functions to map constrained model
            % parameters to unconstrained optimisation variables.
            s0        = linker('s0'       , 'bounded', 0    , 1);
            diffPar   = linker('diffPar'  , 'bounded', 1e-10, 3e-9);
            diffPerp = linker('diffPerp', 'bounded', 1e-10, 3e-9);
            theta     = linker('theta');
            phi       = linker('phi');
            alpha     = linker('alpha');
            kappa1  = linker('kappa1' , 'bounded', 0, 128); 
            kappa2  = linker('kappa2' , 'bounded', 0, 128);
            obj.links = [s0; diffPar; diffPerp;...
                         theta; phi; alpha; kappa1; kappa2];
            obj.links.addConstraint('diffPar>=diffPerp');
            obj.links.addConstraint('kappa2>=kappa1');             
            
            % set hparams and hparams names
            obj.hyperparams = [50; 50];
            obj.hyperparamsName = {'nBinsTheta'; 'nBinsPhi'};
        end%of constructor    
        
        function [sig, out] = synthesize(obj, params, scheme)
            % synthesize(params, scheme, hyperparams) return DW-MR signals.
            %       Input arguments:
            %                params: Numerical column vector of all model
            %                        parameters
            %                scheme: Diffusion scheme
            %
            %      output arguments:
            %                   sig: Numerical column vector of synthesied
            %                        signal
            
            % validate inputs
            validateattributes(params, {'numeric'},...
                {'column', 'nrows', obj.nParams},...
                'cylinderBDA.synthesize', 'params');

            % read params into individual model parameters
            s0          = params(1);   % b0 signal
            diffPar     = params(2);   % diffusivity [s/m^2]
            diffPerp    = params(3);   % diffusivity [s/m^2]
            theta       = params(4);   % elevation angle
            phi         = params(5);   % azimuth angel
            alpha       = params(6);   % angele beween v1_zrot and v2
            kappa1      = params(7);
            kappa2      = params(8);
            
            % read hyper-parameters
            nBinsTheta = obj.hyperparams(1);
            nBinsPhi   = obj.hyperparams(2);
            
            % get Unit Frame
            U = math.getUnitFrame(theta, phi, alpha);
            cB = cylinderBDA.cBingham(-kappa1, -kappa2); % constant for bingmham distribution
            [probs, thisTheta, thisPhi] = cylinderBDA.discritiseBinghamDistribution(U, nBinsTheta, nBinsPhi, kappa1, kappa2, cB);  
            
            % compute cylinder axes for each (theta, phi) [nInstance x 3]
            n = [cos(thisPhi(:)).*sin(thisTheta(:)),...
                 sin(thisPhi(:)).*sin(thisTheta(:)),...
                 cos(thisTheta(:))];
            
            ecDiffPar = sum(((diffPar-diffPerp)*(n*U(:,1)).^2 + diffPerp).*probs.*sin(thisTheta(:)), 1);
            ecDiffPerp1 = sum(((diffPar-diffPerp)*(n*U(:,2)).^2 + diffPerp).*probs.*sin(thisTheta(:)), 1);
            ecDiffPerp2 = diffPar+ 2*diffPerp - ecDiffPar - ecDiffPerp1;

            % gradient matrix [nScheme x 3] for each measurment
            G_dir = [scheme.x, scheme.y, scheme.z];
            
            b = (scheme.DELTA-scheme.delta/3).* ...
                ((scheme.delta .*scheme.G_mag)*math.GAMMA).^2;
            
            s = exp(-b.*(ecDiffPar  *(G_dir*U(:,1)).^2 + ...
                         ecDiffPerp1*(G_dir*U(:,2)).^2 + ...
                         ecDiffPerp2*(G_dir*U(:,3)).^2)); 
            sig = s0 * s;
            if nargout==2
                out.s = s;
                out.ecDiffPar = ecDiffPar;
                out.ecDiffPerp1 = ecDiffPerp1;
                out.ecDiffPerp2 = ecDiffPerp2;
                out.n = n;
                out.U = U;
                out.cB = cB;
                out.probs = probs;
                out.thisTheta = thisTheta;
                out.thisPhi = thisPhi;
                out.b = b;
                out.G_dir = G_dir;
            end
        end
        
         function jac = jacobian(obj, params, scheme)
            % jacobian(params, scheme, hyperparams) return the gradient of 
            % signal wrt to model parameters.
            %
            %       Input arguments:
            %                params: Numerical column vector of all model
            %                        parameters
            %                scheme: Diffusion scheme
            %           hyperparams: non-optimisable model parameters
            %
            %      output arguments:
            %                   jac: Numerical matrix of size M x nParams.
            %
            
            [f0, out] = obj.synthesize(params, scheme);  
            
            % read params into individual model parameters
            s0          = params(1);   % b0 signal
            diffPar     = params(2);   % diffusivity [s/m^2]
            diffPerp    = params(3);   % diffusivity [s/m^2]
            theta       = params(4);   % elevation angle
            phi         = params(5);   % azimuth angel
            alpha       = params(6);   % angele beween v1_zrot and v2
            kappa1      = params(7);
            kappa2      = params(8);
            
            % read hyper-parameters
            nBinsTheta = obj.hyperparams(1);
            nBinsPhi   = obj.hyperparams(2);
            
            % read out parameters
            gf0_gs0 = out.s;
            ecDiffPar = out.ecDiffPar;
            ecDiffPerp1 = out.ecDiffPerp1;
            ecDiffPerp2 = out.ecDiffPerp2;
            n = out.n;
            U = out.U;
            cB = out.cB;
            probs = out.probs;
            thisTheta = out.thisTheta;
            b = out.b;
            G_dir = out.G_dir;
                
            % initialise the jac with zeros
            jac = zeros(size(scheme,1), obj.nParams);
                      
            % gradient wrt to s0
            jac(:,1) = gf0_gs0;
            
            % gradient wrt dPar
            gf_gEcDiffPar = (-b.*((G_dir*U(:,1)).^2)).*f0;
            
            % gradient wrt dPerp1
            gf_gEcDiffPerp1 = (-b.*((G_dir*U(:,2)).^2)).*f0;
            
            % gradient wrt dPerp2
            gf_gEcDiffPerp2 = (-b.*((G_dir*U(:,3)).^2)).*f0;
            
            % gradient of dPar wrt diffPar
            gEcDiffPar_gDiffPar = sum((n*U(:,1)).^2.*probs.*sin(thisTheta(:)), 1);
            gEcDiffPerp1_gDiffPar = sum((n*U(:,2)).^2.*probs.*sin(thisTheta(:)), 1);
            gEcDiffPerp2_gDiffPar = 1 - gEcDiffPar_gDiffPar - gEcDiffPerp1_gDiffPar;
            
            jac(:,2) = gf_gEcDiffPar.*gEcDiffPar_gDiffPar + gf_gEcDiffPerp1.*gEcDiffPerp1_gDiffPar + gf_gEcDiffPerp2.*gEcDiffPerp2_gDiffPar;
            
            % gradient of dPar wrt diffPerp
            gEcDiffPar_gDiffPerp = sum((1-(n*U(:,1)).^2).*probs.*sin(thisTheta(:)), 1);
            gEcDiffPerp1_gDiffPerp = sum((1-(n*U(:,2)).^2).*probs.*sin(thisTheta(:)), 1);
            gEcDiffPerp2_gDiffPerp = 2 - gEcDiffPar_gDiffPerp - gEcDiffPerp1_gDiffPerp;
            
            jac(:,3) = gf_gEcDiffPar.*gEcDiffPar_gDiffPerp + gf_gEcDiffPerp1.*gEcDiffPerp1_gDiffPerp + gf_gEcDiffPerp2.*gEcDiffPerp2_gDiffPerp;
            
            % gradient wrt theta, phi, and alpha
            % using chain rule, comput gf_gV1, gf_gV2, gf_gV3.           
            gDPar_gU1 = sum(2*(diffPar-diffPerp)*(n*U(:,1)).*n.*probs.*sin(thisTheta(:)), 1);
            gDPar_gU2 = zeros(1,3);
            gDPar_gU3 = zeros(1,3);
            
            gDPerp1_gU1 = zeros(1,3);
            gP_gU2 = (-2*kappa1)*probs.*(n*U(:,2)).*n;
            gP_gU3 = (-2*kappa2)*probs.*(n*U(:,3)).*n;
            gDPerp1_gU2 = sum(2*(diffPar-diffPerp)*(n*U(:,2)).*n.*probs.*sin(thisTheta(:)), 1) + ...
                          sum(((diffPar-diffPerp)*(n*U(:,2)).^2 + diffPerp).*sin(thisTheta(:)).*gP_gU2, 1);
            gDPerp1_gU3 = sum(((diffPar-diffPerp)*(n*U(:,2)).^2 + diffPerp).*sin(thisTheta(:)).*gP_gU3, 1);
            
            gDPerp2_gU1 = - gDPar_gU1 - gDPerp1_gU1;
            gDPerp2_gU2 = - gDPar_gU2 - gDPerp1_gU2;
            gDPerp2_gU3 = - gDPar_gU3 - gDPerp1_gU3;
            
            gf_gU1 = repmat(-2*ecDiffPar*b.*(G_dir*U(:,1)).*f0,1,3).*G_dir   + ...
                     repmat(gf_gEcDiffPar, 1, 3).*gDPar_gU1          + ...
                     repmat(gf_gEcDiffPerp1, 1, 3).*gDPerp1_gU1      + ...
                     repmat(gf_gEcDiffPerp2, 1, 3).*gDPerp2_gU1;
                 
            gf_gU2 = repmat(-2*ecDiffPerp1*b.*(G_dir*U(:,2)).*f0,1,3).*G_dir + ...
                     repmat(gf_gEcDiffPar, 1, 3).*gDPar_gU2          + ...
                     repmat(gf_gEcDiffPerp1, 1, 3).*gDPerp1_gU2      + ...
                     repmat(gf_gEcDiffPerp2, 1, 3).*gDPerp2_gU2;
                 
            gf_gU3 = repmat(-2*ecDiffPerp2*b.*(G_dir*U(:,3)).*f0,1,3).*G_dir + ...
                     repmat(gf_gEcDiffPar, 1, 3).*gDPar_gU3          + ...
                     repmat(gf_gEcDiffPerp1, 1, 3).*gDPerp1_gU3      + ...
                     repmat(gf_gEcDiffPerp2, 1, 3).*gDPerp2_gU3;
            
            % gradient wrt theta
            gU1_gTheta = [ cos(phi)*cos(theta); ...
                           sin(phi)*cos(theta); ...
                          -sin(theta)];
                      
            gU2_gTheta = [-cos(alpha)*cos(phi)*sin(theta); ...
                          -cos(alpha)*sin(phi)*sin(theta); ...
                          -cos(alpha)*cos(theta)];
                     
            gU3_gTheta = [ sin(theta)*cos(phi)*sin(alpha);...
                           sin(theta)*sin(phi)*sin(alpha);...
                           cos(theta)*sin(alpha)];
                       
            jac(:,4) = gf_gU1*gU1_gTheta + gf_gU2*gU2_gTheta + gf_gU3*gU3_gTheta;          
            
            % gradient wrt phi
            gV1_gPhi   = [-sin(phi)*sin(theta); ...
                           cos(phi)*sin(theta); ...
                           0];
                       
            gV2_gPhi   = [-sin(phi)*cos(alpha)*cos(theta) - cos(phi)*sin(alpha); ...          
                           cos(phi)*cos(alpha)*cos(theta) - sin(phi)*sin(alpha); ...
                           0];
                       
            gV3_gPhi   = [ sin(phi)*cos(theta)*sin(alpha) - cos(phi)*cos(alpha);...
                          -cos(phi)*cos(theta)*sin(alpha) - sin(phi)*cos(alpha);...
                         0];    
                     
            jac(:,5) = gf_gU1*gV1_gPhi   + gf_gU2*gV2_gPhi   + gf_gU3*gV3_gPhi;
            
            % gradient wrt alpha
            % gV1_gAlpha = zeros(3, 1);             
            
            gV2_gAlpha = [-sin(alpha)*cos(theta)*cos(phi) - cos(alpha)*sin(phi);...
                          -sin(alpha)*cos(theta)*sin(phi) + cos(alpha)*cos(phi);
                           sin(alpha)*sin(theta)];                     
           
            gV3_gAlpha = [-cos(alpha)*cos(theta)*cos(phi)+sin(alpha)*sin(phi);...
                          -cos(alpha)*cos(theta)*sin(phi)-sin(alpha)*cos(phi);...
                           cos(alpha)*sin(theta)];
           
            jac(:,6) = gf_gU2*gV2_gAlpha + gf_gU3*gV3_gAlpha;
            
            [gcB_gK1, gcB_gK2] = cylinderBDA.cBinghamGrad(-kappa1, -kappa2);
            
            gP_gK1      = probs.*(gcB_gK1/cB-(n*U(:,2)).^2);
            gP_gK2      = probs.*(gcB_gK2/cB-(n*U(:,3)).^2);
            
            gDPar_gK1 = sum(((diffPar-diffPerp)*(n*U(:,1)).^2+ diffPerp).*gP_gK1.*sin(thisTheta(:)), 1);
            gDPerp1_gK1 = sum(((diffPar-diffPerp)*(n*U(:,2)).^2+ diffPerp).*gP_gK1.*sin(thisTheta(:)), 1);
            gDPerp2_gK1 = - gDPar_gK1 - gDPerp1_gK1;
            
            jac(:,7) = gf_gEcDiffPar.*gDPar_gK1 + gf_gEcDiffPerp1.*gDPerp1_gK1 + ...
                       gf_gEcDiffPerp2.*gDPerp2_gK1;
            
            gDPar_gK2 = sum(((diffPar-diffPerp)*(n*U(:,1)).^2+ diffPerp).*gP_gK2.*sin(thisTheta(:)), 1);
            gDPerp1_gK2 = sum(((diffPar-diffPerp)*(n*U(:,2)).^2+ diffPerp).*gP_gK2.*sin(thisTheta(:)), 1);
            gDPerp2_gK2 = - gDPar_gK2 - gDPerp1_gK2;
            
            jac(:,8) = gf_gEcDiffPar.*gDPar_gK2 + gf_gEcDiffPerp1.*gDPerp1_gK2 + ...
                       gf_gEcDiffPerp2.*gDPerp2_gK2;       
        end % of jacobian
        %\\
    end % of method (public) 
   %
end%of class zeppelinBDA
