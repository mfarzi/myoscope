classdef cylinderBDA < compartment
    % CYLINDERBDA
    % (Cylinder Model with Bingham Distribution of Axes)
    %
    %   A CYLINDERBDA object is a basic COMPARTMENT object represenitng 
    %   particles diffusion in cylinders with various longitudinal axis 
    %   drawn from a Bingham distribution. This class implements inherited
    %   abstract methods and properties. This class should be used as an 
    %   input to the MULTICOMPARTMENT class.
    %
    %   Model Parameters:
    %       s0                - Normalised b0-signal 
    %       diffPar           - diffusivity along the cylinder axis [m^2/s]
    %       r                 - radius of cylinder [m]
    %       theta             - elevation angle; the angle between the
    %                           cylinder axis (n1) and the z-axis. [radian]
    %       phi               - azimuth angle; the angle between the x-axis
    %                           and the n1 projection onto the xy-plane. 
    %                           [radian]
    %       alpha             - plane agnel; the angle between the second
    %                           eigenvector and the v1 rotated by pi/2 
    %                           around the z-axis. [radian]
    %       kappa1            - Concentration paramter for  the Bingham  
    %                           distribution
    %       kappa2            - Concentration paramter for the Bingham 
    %                           distribution                            
    %
    %   For mathematical background see
    %       Y. Assaf, R. Z. Freidlin, G. K. Rohde, and P. J. Basser,
    %       New modeling and experimental framework to characterize 
    %       hindered and restricted water diffusion in brain white 
    %       matter, Magn. Reson. Med., vol. 52, no. 5, pp. 965?978, 2004.
    %
    %   See also: multicompartment, cylinder, cylinderECS, cylinderGDR, 
    %             stick
    %
    % Mohsen Farzi
    % Email: m.farzi@leeds.ac.uk
    
    properties (SetAccess='protected')
        name = 'cylinderBDA';    
        nCompartments = 1;     % number of basic COMPARTMETNT objects 
        nParams = 8;
        nHyperparams = 2;
    end
    
    properties (Access = 'protected')
        comp = [];             % List of COMPARTMENT object 
        links;
        hyperparams = [];      % non-optimisable model parameters
        hyperparamsName = {};
    end
    
    methods 
        function obj = cylinderBDA()
            %CYLINDERBDA Construct Function.
            %
            %   cylinderBDA() constructs a single-compartment model
            
            % initialise the linking functions to map constrained model
            % parameters to unconstrained optimisation variables.
            s0      = linker('s0'     , 'bounded', 0    , 1);            
            diffPar = linker('diffPar', 'bounded', 1e-10, 3e-9);               
            r       = linker('r'      , 'bounded', 1e-6 , 20e-6);       
            theta   = linker('theta'); 
            phi     = linker('phi');
            alpha   = linker('alpha');
            kappa1  = linker('kappa1' , 'bounded', 0, 128); 
            kappa2  = linker('kappa2' , 'bounded', 0, 128);
            
            obj.links = [s0; diffPar; r; theta; phi; alpha; kappa1; kappa2];
            obj.links.addConstraint('kappa2>=kappa1');
            
            % set hparams and hparams names
            obj.hyperparams = [50; 50];
            obj.hyperparamsName = {'nBinsTheta'; 'nBinsPhi'};
        end
        
        function [s, tmpAccessMemory] = synthesize(obj, params, scheme)
            % synthesize(params, scheme, hyperparams) return DW-MR signals.
            %       Input arguments:
            %                params: Numerical column vector of all model
            %                        parameters
            %                scheme: Diffusion scheme
            %           hyperparams: Integer column vector [2x1]
            %
            %      output arguments:
            %                     s: Numerical column vector of synthesied
            %                        signal
            
            % validate inputs
            validateattributes(params, {'numeric'},...
                {'column', 'nrows', obj.nParams},...
                'cylinderBDA.synthesize', 'params');

            % read params into individual model parameters
            s0          = params(1);   % b0 signal
            diffPar     = params(2);   % diffusivity [s/m^2]
            r           = params(3);   % cylinder radius [m]
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
            
            % compute cylinder axes for each (theta, phi) [3 x nInstance]
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
            Lpar = cylinder.get_Lpar(scheme, diffPar);
            sPar = exp(-Lpar.*Gpar_mag2);
            
            % compute the signal perpendicular to the cylinder axis
            % [1 x length(Jp1ROOTS)]
            [Lperp, nom, denom] = cylinder.get_Lperp(scheme, diffPar, r);
            sPerp = exp(-Lperp.*Gperp_mag2);
            
            % compute signal as the product of the value along and
            % perpendicular to the cylinder axis
            sig = s0 * sPerp .* sPar;
            
            % compute the average wrt Bingham Distribution
            weights = probs.*sin(thisTheta(:));
            s = sig*weights;
            gf0_gs0 = (sPerp.*sPar)*weights;
            
            if nargout == 2
                tmpAccessMemory.gf0_gs0 = gf0_gs0;
                tmpAccessMemory.sig   = sig;
                tmpAccessMemory.Gpar_mag2 = Gpar_mag2;
                tmpAccessMemory.Gperp_mag2 = Gperp_mag2;
                tmpAccessMemory.Lpar  = Lpar;
                tmpAccessMemory.Lperp = Lperp;
                tmpAccessMemory.nom   = nom;
                tmpAccessMemory.denom = denom;
                tmpAccessMemory.probs = probs;
                tmpAccessMemory.weights = weights;
                tmpAccessMemory.theta = thisTheta(:);
                tmpAccessMemory.phi   = thisPhi(:);
                tmpAccessMemory.n     = n;
                tmpAccessMemory.U = U;
                tmpAccessMemory.cB = cB;
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
            %
            %      output arguments:
            %                   jac: Numerical matrix of size M x nParams.
            %
            % synthesize signal
            [f0, tmpAccessMemory] = obj.synthesize(params, scheme);
            
            % read intermediate variables
            gf0_gs0 = tmpAccessMemory.gf0_gs0;
            nom = tmpAccessMemory.nom;
            denom = tmpAccessMemory.denom;
            Gpar_mag2 = tmpAccessMemory.Gpar_mag2;
            Gperp_mag2 = tmpAccessMemory.Gperp_mag2;
            sig = tmpAccessMemory.sig;
            probs = tmpAccessMemory.probs;
            weights = tmpAccessMemory.weights;
            thisTheta = tmpAccessMemory.theta;
            thisPhi = tmpAccessMemory.phi;
            U = tmpAccessMemory.U;
            cB = tmpAccessMemory.cB;
            
            % read params into individual model parameters
            s0          = params(1);   % b0 signal
            diffPar     = params(2);   % diffusivity [s/m^2]
            r           = params(3);   % cylinder radius [m]
            theta       = params(4);   % elevation angle
            phi         = params(5);   % azimuth angel
            alpha       = params(6);   % angele beween v1_zrot and v2
            kappa1      = params(7);
            kappa2      = params(8);
            
            % initialise the jac with zeros
            jac = zeros(size(scheme, 1), obj.nParams);
            
            %%
            % gSgS0 : gradient wrt s0
            jac(:,1) = gf0_gs0;
            
            % gradient wrt diffPar
            % first compute gradients for each pair f(theta', phi') and then
            % average over the Bingham distriubtion.
            gLpar_gDiffPar = cylinder.get_gLpar_gDiff(scheme);
            gLperp_gDiffPar = cylinder.get_gLperp_gDiff(scheme, diffPar, r, nom, denom);
            
            % gradient wrt diffPar      
            gf_gDiffPar = -(gLpar_gDiffPar.*Gpar_mag2 + gLperp_gDiffPar.*Gperp_mag2).*sig;
            
            % gSgDiffPar: average gradient over the Bingham distribution
            jac(:,2) = gf_gDiffPar*weights; 
            
            % gradient wrt r
            % first compute gradients for each pair f(theta', phi') and then
            % average over the Bingham distriubtion.
            gLperp_gR = cylinder.get_gLperp_gR(scheme, diffPar, r, nom, denom);
            gf_gR = -(gLperp_gR.*Gperp_mag2).*sig;
            
            % gSgR
            jac(:,3)  = gf_gR*weights;
            
            %%          
            u2 = U(:,2); u3 = U(:,3);
            
            [nX, nY, nZ]  = cylinderBDA.sphericalToCanonical(thisTheta, thisPhi, 1);
            u2DotN = u2(1)*nX + u2(2)*nY + u2(3)*nZ;
            u3DotN = u3(1)*nX + u3(2)*nY + u3(3)*nZ;
            
            gFb_gU2 = -2*kappa1*u2DotN.*probs.*[nX, nY, nZ];
            gFb_gU3 = -2*kappa2*u3DotN.*probs.*[nX, nY, nZ];
            
            % gradient wrt theta
            gU2_gTheta = [-cos(alpha)*cos(phi)*sin(theta); ...
                          -cos(alpha)*sin(phi)*sin(theta); ...
                          -cos(alpha)*cos(theta)];
                     
            gU3_gTheta = [ sin(theta)*cos(phi)*sin(alpha);...
                           sin(theta)*sin(phi)*sin(alpha);...
                           cos(theta)*sin(alpha)];
                       
            gFb_gTheta = gFb_gU2*gU2_gTheta + gFb_gU3*gU3_gTheta;   
            % gSgTheta 
            jac(:,4)   = sig*(gFb_gTheta.*sin(thisTheta));
            
            % gradient wrt phi
            gU2_gPhi   = [-sin(phi)*cos(alpha)*cos(theta) - cos(phi)*sin(alpha); ...          
                           cos(phi)*cos(alpha)*cos(theta) - sin(phi)*sin(alpha); ...
                           0];
                       
            gU3_gPhi   = [ sin(phi)*cos(theta)*sin(alpha) - cos(phi)*cos(alpha);...
                          -cos(phi)*cos(theta)*sin(alpha) - sin(phi)*cos(alpha);...
                         0];    
                     
            gFb_gPhi = gFb_gU2*gU2_gPhi + gFb_gU3*gU3_gPhi;
            % gSgPhi
            jac(:,5)   = sig*(gFb_gPhi.*sin(thisTheta));
            
            % gradient wrt alpha             
            gU2_gAlpha = [-sin(alpha)*cos(theta)*cos(phi) - cos(alpha)*sin(phi);...
                          -sin(alpha)*cos(theta)*sin(phi) + cos(alpha)*cos(phi);
                           sin(alpha)*sin(theta)];                     
           
            gU3_gAlpha = [-cos(alpha)*cos(theta)*cos(phi)+sin(alpha)*sin(phi);...
                          -cos(alpha)*cos(theta)*sin(phi)-sin(alpha)*cos(phi);...
                           cos(alpha)*sin(theta)];
                       
            gFb_gAlpha = gFb_gU2*gU2_gAlpha + gFb_gU3*gU3_gAlpha;           
            % gSgAlpha
            jac(:,6)   = sig*(gFb_gAlpha.*sin(thisTheta));
            
            [gcB_gK1, gcB_gK2] = cylinderBDA.cBinghamGrad(-kappa1, -kappa2);
            % gSgKappa1
            jac(:,7)   = sig*(weights.*(-u2DotN.^2+gcB_gK1/cB));
            
            % gSgKappa2
            jac(:,8)   = sig*(weights.*(-u3DotN.^2+gcB_gK2/cB));
        end
        
    end%of methods (public)
        
    methods (Static)
        function [probs, thisTheta, thisPhi] = discritiseBinghamDistribution(U, nBinsTheta, nBinsPhi, kappa1, kappa2, cB)
            thetaVec = linspace(0, pi, nBinsTheta+1);
            phiVec = linspace(0, 2*pi, nBinsPhi+1);
            
            [thisTheta, thisPhi] = meshgrid((thetaVec(1:nBinsTheta)+thetaVec(2:nBinsTheta+1))/2,...
                (phiVec(1:nBinsPhi)+phiVec(2:nBinsPhi+1))/2);
            
            dS = (thetaVec(2:nBinsTheta+1)-thetaVec(1:nBinsTheta))'*...
                 (phiVec(2:nBinsPhi+1)-phiVec(1:nBinsPhi));
            dS = abs(dS/(4*pi));
            
            [x, y, z] = cylinderBDA.sphericalToCanonical(thisTheta, thisPhi, 1);
            
            n = [x(:), y(:), z(:)];
            
            probs = cylinderBDA.binghamPdf(n', U(:,2), U(:,3), -kappa1, -kappa2, cB);
            probs = probs(:).*dS(:);
        end
        
        function p = binghamPdf(n, u1, u2, kappa1, kappa2, cB)
            p = exp(kappa1*(u1'*n).^2 + kappa2*(u2'*n).^2)/cB;
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
    end%of methods (private)
    %\\
end%of cylinderBDA class

