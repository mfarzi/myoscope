classdef cylinderECS < compartment
    % CYLINDERECS 
    % (Cylinder with Elliptcial Cross Section)
    %
    %   A CYLINDERECS object is a basic COMPARTMENT object represenitng
    %   particles diffusion in an elliptical cylinder. This class 
    %   implements inherited abstract methods and properties. This class 
    %   should be used as input to MULTICOMPARTMENT class.
    %
    %   Model Parameters:
    %       s0                - Normalised b0-signal 
    %       diffPar           - diffusivity along the cylinder axis [m^2/s]
    %       r1                - radius of ellipsoid along the major axis
    %                           [m]
    %       r2                - radius of ellipsoid along the minor axis
    %                           [m]
    %       theta             - elevation angle; the angle between the
    %                           cylinder axis (n1) and the z-axis. [radian]
    %       phi               - azimuth angle; the angle between the x-axis
    %                           and the n1 projection onto the xy-plane. 
    %                           [radian]
    %       alpha             - the angle between the major ellipsoid axis
    %                           and the n1 projection onto the cross 
    %                           sectional elliposid plain [radian]
    %
    %   For mathematical background see
    %       Vangelderen, P., DesPres, D., Vanzijl, P.C.M. and Moonen, C.,
    %       "Evaluation of restricted diffusion in cylinders.
    %       Phosphocreatine in rabbit leg muscle.", Journal of Magnetic
    %       Resonance, Series B, 103(3), pp.255-260, 1994.
    %
    %   See also: multicompartment, cylinder, cylinderGDR, cylinderBDA,
    %             stick
    %
    % Mohsen Farzi
    % Email: m.farzi@leeds.ac.uk
    
    properties (SetAccess='protected')
        name = 'cylinderECS';    
        nCompartments = 1;     % number of basic COMPARTMETNT objects 
        nParams = 7;
        nHyperparams = 0;
    end
    
    properties (Access = 'protected')
        comp = [];             % List of COMPARTMENT object 
        links;
        hyperparams = [];
        hyperparamsName = {};
    end   
    
    methods
        function obj = cylinderECS()
            %cylinderECS Construct Function.
            %
            %   cylinderECS() constructs a single-compartment model
            
            % initialise the linking functions to map constrained model
            % parameters to unconstrained optimisation variables.
            s0        = linker('s0'       , 'bounded', 0    , 1);
            diffPar   = linker('diffPar'  , 'bounded', 1e-10, 3e-9);
            r1 = linker('r1', 'bounded', 1e-6, 20e-6);
            r2 = linker('r2', 'bounded', 1e-6, 20e-6);
            theta     = linker('theta');
            phi       = linker('phi');
            alpha     = linker('alpha');
            obj.links = [s0; diffPar; r1; r2; theta; phi; alpha];
            obj.links.addConstraint('r1>=r2');
        end
        
        function [sig, out] = synthesize(obj, params, scheme)
            % synthesize(params, scheme, hyperparams) return DW-MR signals
            % from a cylinderECS. The signal is modeled as the product of 
            % signals parallel and perpendicular to the cylinder axis.
            %       Input arguments:
            %                params: Numerical column vector of all model
            %                        parameters
            %                scheme: Diffusion scheme
            %
            %      output arguments:
            %                     s: Numerical column vector of synthesied
            %                        signal
            %  
            
            % validate inputs
            validateattributes(params, {'double'},...
                {'column', 'nrows', obj.nParams},...
                'cylinderECS.synthesize', 'params');
            
             % read params into individual model parameters
            s0      = params(1);   % b0 signal
            diffPar = params(2);   % diffusivity [s/m^2]
            r1      = params(3);   % radius along the long axis [m]
            r2      = params(4);   % radius along the short axis [m]
            theta   = params(5);   % elevation angle [radian]
            phi     = params(6);   % azimuth angel [radian]
            alpha   = params(7);   % plane angle [radian]
            
            
            G = [scheme.x, scheme.y, scheme.z].*scheme.G_mag;
            N = math.getUnitFrame(theta, phi, alpha); 
            n1 = N(:,1); n2 = N(:,2); n3 = N(:,3);
            
            % comput the gradient along and perpendicular
            % to the cylinder axis 
            Gpar_mag2 = (G*n1).^2; 
            Gperp1_mag2 = (G*n2).^2;
            Gperp2_mag2 = (G*n3).^2;
            
            % compute bPar
            bPar = (scheme.DELTA-scheme.delta/3).*(math.GAMMA*scheme.delta).^2*diffPar;
            
            % compute bPerp1
            beta_r1 = math.Jp1ROOTS/r1;          
            DB2_r1 = diffPar*beta_r1.^2;
            nom1 = 2*scheme.delta*DB2_r1' - 2  ...
                + 2*exp(-scheme.delta*DB2_r1') ...
                + 2*exp(-scheme.DELTA*DB2_r1') ...
                - exp(-(scheme.DELTA-scheme.delta)*DB2_r1') ...
                - exp(-(scheme.DELTA+scheme.delta)*DB2_r1');

            denom1 = diffPar^2*beta_r1.^6.*(math.Jp1ROOTS.^2-1);
            
            bPerp1 = 2*math.GAMMA^2*nom1*(1./denom1);
            
            % compute bPerp2
            beta_r2 = math.Jp1ROOTS/r2;          
            DB2_r2 = diffPar*beta_r2.^2;
            nom2 = 2*scheme.delta*DB2_r2' - 2  ...
                + 2*exp(-scheme.delta*DB2_r2') ...
                + 2*exp(-scheme.DELTA*DB2_r2') ...
                - exp(-(scheme.DELTA-scheme.delta)*DB2_r2') ...
                - exp(-(scheme.DELTA+scheme.delta)*DB2_r2');

            denom2 = diffPar^2*beta_r2.^6.*(math.Jp1ROOTS.^2-1);
            
            bPerp2 = 2*math.GAMMA^2*nom2*(1./denom2);
            
            
            % compute signal as the product of the value along and
            % perpendicular to the cylinder axis
            s = exp(-bPar  .*Gpar_mag2)   .* ...
                exp(-bPerp1.*Gperp1_mag2) .* ...
                exp(-bPerp2.*Gperp2_mag2);
            sig = s0*s;
            if nargout==2
                out.s = s;
                out.n1 = n1;
                out.n2 = n2;
                out.n3 = n3;
                out.G = G;
                out.Gpar_mag2 = Gpar_mag2;
                out.Gperp1_mag2 = Gperp1_mag2;
                out.Gperp2_mag2 = Gperp2_mag2;
                out.bPar = bPar;
                out.bPerp1 = bPerp1;
                out.bPerp2 = bPerp2;
                out.nom1 = nom1;
                out.denom1 = denom1;
                out.nom2 = nom2;
                out.denom2 = denom2;
            end
        end % of synthesize
        
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
            %                   jac: Numerical matrix of size Mx7.
            % 
            
            [f0, out] = obj.synthesize(params, scheme, []);
            
            % read params into individual model parameters
            s0      = params(1);   % b0 signal
            diffPar = params(2);   % diffusivity [s/m^2]
            r1      = params(3);   % radius along the long axis [m]
            r2      = params(4);   % radius along the short axis [m]
            theta   = params(5);   % elevation angle [radian]
            phi     = params(6);   % azimuth angel [radian]
            alpha   = params(7);   % plane angle [radian]
            
            % read intermediate variables
            gf0_gs0 = out.s;
            n1 = out.n1;
            n2 = out.n2;
            n3 = out.n3;
            G = out.G;
            Gpar_mag2 = out.Gpar_mag2;
            Gperp1_mag2 = out.Gperp1_mag2;
            Gperp2_mag2 = out.Gperp2_mag2;
            bPar = out.bPar;
            bPerp1 = out.bPerp1;
            bPerp2 = out.bPerp2;
            nom1 = out.nom1;
            denom1 = out.denom1;
            nom2 = out.nom2;
            denom2 = out.denom2;
            
            
            % initialise the jac with zeros
            jac = zeros(size(scheme, 1), obj.nParams);
           
            
            % compute bPerp1
            beta_r1 = math.Jp1ROOTS/r1;          
            B2_r1 = beta_r1.^2;
            
            % compute bPerp2
            beta_r2 = math.Jp1ROOTS/r2;          
            B2_r2 = beta_r2.^2;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%       Gradient wrt to s0-jac(:,1)        %%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % gradient wrt to s0
            jac(:,1) = gf0_gs0;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%       Gradient wrt to diffPar-jac(:,2)      %%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % gbPar_gDiffPar
            gbPar_gDiffPar = bPar/diffPar;
            
            % gbPerp1_gDiffPar
            gNom1_gDiffPar = 2*scheme.delta*B2_r1' ...
                        - 2*exp(-diffPar*scheme.delta*B2_r1').*(scheme.delta*B2_r1') ...
                        - 2*exp(-diffPar*scheme.DELTA*B2_r1').*(scheme.DELTA*B2_r1') ...
                        +   exp(-diffPar*(scheme.DELTA-scheme.delta)*B2_r1').*((scheme.DELTA-scheme.delta)*B2_r1')...
                        +   exp(-diffPar*(scheme.DELTA+scheme.delta)*B2_r1').*((scheme.DELTA+scheme.delta)*B2_r1');
            
            gDenom1_gDiffPar = 2*diffPar*beta_r1.^6.*(math.Jp1ROOTS.^2-1);
            tmp = bsxfun(@times, gNom1_gDiffPar, denom1') - bsxfun(@times, nom1, gDenom1_gDiffPar');
            tmp = bsxfun(@rdivide, tmp, denom1'.^2);
            gbPerp1_gDiffPar = sum(tmp, 2)*2*math.GAMMA^2;
     
            % gbPerp2_gDiffPar
            gNom2_gDiffPar = 2*scheme.delta*B2_r2' ...
                        - 2*exp(-diffPar*scheme.delta*B2_r2').*(scheme.delta*B2_r2') ...
                        - 2*exp(-diffPar*scheme.DELTA*B2_r2').*(scheme.DELTA*B2_r2') ...
                        +   exp(-diffPar*(scheme.DELTA-scheme.delta)*B2_r2').*((scheme.DELTA-scheme.delta)*B2_r2')...
                        +   exp(-diffPar*(scheme.DELTA+scheme.delta)*B2_r2').*((scheme.DELTA+scheme.delta)*B2_r2');
            
            gDenom2_gDiffPar = 2*diffPar*beta_r2.^6.*(math.Jp1ROOTS.^2-1);
            tmp = bsxfun(@times, gNom2_gDiffPar, denom2') - bsxfun(@times, nom2, gDenom2_gDiffPar');
            tmp = bsxfun(@rdivide, tmp, denom2'.^2);
            gbPerp2_gDiffPar = sum(tmp, 2)*2*math.GAMMA^2;
            
            % gradient wrt diffPar       
            gf_gDiffPar = -(gbPar_gDiffPar.*Gpar_mag2     + ...
                         gbPerp1_gDiffPar.*Gperp1_mag2 + ...
                         gbPerp2_gDiffPar.*Gperp2_mag2).*f0;
            jac(:,2) = gf_gDiffPar;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%       Gradient wrt to r1-jac(:,3)        %%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % gbPar_gR1 = 0;
            
            % gbPerp1_gR1
            gNom1_gBm = 4*diffPar*scheme.delta*beta_r1' ...
                      - 2*exp(-diffPar*scheme.delta*B2_r1').*(2*diffPar*scheme.delta*beta_r1') ...
                      - 2*exp(-diffPar*scheme.DELTA*B2_r1').*(2*diffPar*scheme.DELTA*beta_r1') ...
                      +   exp(-diffPar*(scheme.DELTA-scheme.delta)*B2_r1').*(2*diffPar*(scheme.DELTA-scheme.delta)*beta_r1') ...
                      +   exp(-diffPar*(scheme.DELTA+scheme.delta)*B2_r1').*(2*diffPar*(scheme.DELTA+scheme.delta)*beta_r1');
            
            gDenom1_gBm = 6*diffPar^2*beta_r1.^5.*(math.Jp1ROOTS.^2-1);
            
            tmp = bsxfun(@times, gNom1_gBm, denom1') - bsxfun(@times, nom1, gDenom1_gBm');
            gbPerp1_gBm = 2*math.GAMMA^2*bsxfun(@rdivide, tmp, denom1'.^2);
            
            % gbPerp1_gR1
            gBm_gR1 = -math.Jp1ROOTS/r1^2;
            gbPerp1_gR1 = gbPerp1_gBm * gBm_gR1;
            
            % gradient wrt R1
            gf_gR1 = -(gbPerp1_gR1.*Gperp1_mag2).*f0;
            jac(:,3) = gf_gR1;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%       Gradient wrt to r2-jac(:,4)        %%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % gbPerp2_gR
            gNom2_gBm = 4*diffPar*scheme.delta*beta_r2' ...
                      - 2*exp(-diffPar*scheme.delta*B2_r2').*(2*diffPar*scheme.delta*beta_r2') ...
                      - 2*exp(-diffPar*scheme.DELTA*B2_r2').*(2*diffPar*scheme.DELTA*beta_r2') ...
                      +   exp(-diffPar*(scheme.DELTA-scheme.delta)*B2_r2').*(2*diffPar*(scheme.DELTA-scheme.delta)*beta_r2') ...
                      +   exp(-diffPar*(scheme.DELTA+scheme.delta)*B2_r2').*(2*diffPar*(scheme.DELTA+scheme.delta)*beta_r2');
            
            gDenom2_gBm = 6*diffPar^2*beta_r2.^5.*(math.Jp1ROOTS.^2-1);
            
            tmp = bsxfun(@times, gNom2_gBm, denom2') - bsxfun(@times, nom2, gDenom2_gBm');
            gbPerp2_gBm = 2*math.GAMMA^2*bsxfun(@rdivide, tmp, denom2'.^2);
            
            % gbPerp2_gR2
            gBm_gR2 = -math.Jp1ROOTS/r2^2;
            gbPerp2_gR2 = gbPerp2_gBm * gBm_gR2;
            
            % gradient wrt R2
            gf_gR2 = -(gbPerp2_gR2.*Gperp2_mag2).*f0;
            jac(:,4) = gf_gR2;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%      Gradient wrt to theta-jac(:,5)      %%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            gf_gN1 = repmat(-2*bPar.*(G*n1).*f0,1,3).*G;
            gf_gN2 = repmat(-2*bPerp1.*(G*n2).*f0,1,3).*G;
            gf_gN3 = repmat(-2*bPerp2.*(G*n3).*f0,1,3).*G;
            
            gN1_gTheta = [ cos(phi)*cos(theta); ...
                           sin(phi)*cos(theta); ...
                          -sin(theta)];
                      
            gN2_gTheta = [-cos(alpha)*cos(phi)*sin(theta); ...
                          -cos(alpha)*sin(phi)*sin(theta); ...
                          -cos(alpha)*cos(theta)];
                     
            gN3_gTheta = [ sin(theta)*cos(phi)*sin(alpha);...
                           sin(theta)*sin(phi)*sin(alpha);...
                           cos(theta)*sin(alpha)];
                       
            jac(:,5) = gf_gN1*gN1_gTheta + gf_gN2*gN2_gTheta + gf_gN3*gN3_gTheta;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%      Gradient wrt to phi-jac(:,6)         %%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            gN1_gPhi   = [-sin(phi)*sin(theta); ...
                           cos(phi)*sin(theta); ...
                           0];
                       
            gN2_gPhi   = [-sin(phi)*cos(alpha)*cos(theta) - cos(phi)*sin(alpha); ...          
                           cos(phi)*cos(alpha)*cos(theta) - sin(phi)*sin(alpha); ...
                           0];
                       
            gN3_gPhi   = [ sin(phi)*cos(theta)*sin(alpha) - cos(phi)*cos(alpha);...
                          -cos(phi)*cos(theta)*sin(alpha) - sin(phi)*cos(alpha);...
                         0];    
                     
            jac(:,6) = gf_gN1*gN1_gPhi   + gf_gN2*gN2_gPhi   + gf_gN3*gN3_gPhi;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%      Gradient wrt to alpha-jac(:,7)       %%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            gN2_gAlpha = [-sin(alpha)*cos(theta)*cos(phi) - cos(alpha)*sin(phi);...
                          -sin(alpha)*cos(theta)*sin(phi) + cos(alpha)*cos(phi);
                           sin(alpha)*sin(theta)];                     
           
            gN3_gAlpha = [-cos(alpha)*cos(theta)*cos(phi)+sin(alpha)*sin(phi);...
                          -cos(alpha)*cos(theta)*sin(phi)-sin(alpha)*cos(phi);...
                           cos(alpha)*sin(theta)];
           
            jac(:,7) = gf_gN2*gN2_gAlpha + gf_gN3*gN3_gAlpha;                 
        end % of jacobian
    end% of methods (public)
    %\\
end % of CylinderECS class
