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
    %   M. Farzi, D. McClymont, H. Whittington, M.C. Zdora, L. Khazin,
    %   C.A. Lygate, C. Rau,E. Dall?Armellina, I. Teh, and J.E. Schneider,
    %   "Assessing myocardial microstructure withbiophysical models of
    %   diffusion MRI," in IEEE Transactions on Medical Imaging,
    %   July 2021, DOI: 10.1109/TMI.2021.3097907 
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
        
        function [sig, out] = synthesize(obj, params, schemefile)
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
            
            assert(isa(schemefile, 'scheme'), ...
                'MATLAB:tensor:invalidInputArgument',...
                'Scheme file should be of type scheme.');
            
            assert(strcmp(schemefile.type, 'stejskal-tanner'), ...
                'MATLAB:tensor:invalidInputArgument',...
                'Scheme file should be of type stejskal-tanner.');
            
            % read params into individual model parameters
            s0      = params(1);   % b0 signal
            diffPar = params(2);   % diffusivity [s/m^2]
            r1      = params(3);   % radius along the long axis [m]
            r2      = params(4);   % radius along the short axis [m]
            theta   = params(5);   % elevation angle [radian]
            phi     = params(6);   % azimuth angel [radian]
            alpha   = params(7);   % plane angle [radian]
            
            
            G = schemefile.ghat.*schemefile.gmag;
            U = math.getUnitFrame(theta, phi, alpha); 
            
            % comput the gradient along and perpendicular
            % to the cylinder axis 
            Gpar_mag2 = (G*U(:,1)).^2; 
            Gperp1_mag2 = (G*U(:,2)).^2;
            Gperp2_mag2 = (G*U(:,3)).^2;
            
            % compute the signal parallel to the cylinder axis
            Lpar = cylinder.get_Lpar(schemefile, diffPar);
            Spar = exp(-Lpar.*Gpar_mag2);
            
            
            % compute the signal perpendicular to the cylinder axis
            % major axis
            [Lperp1, nom1, denom1] = cylinder.get_Lperp(schemefile, diffPar, r1);
            Sperp1 = exp(-Lperp1.*Gperp1_mag2);
            
            % compute the signal perpendicular to the cylinder axis
            % major axis
            [Lperp2, nom2, denom2] = cylinder.get_Lperp(schemefile, diffPar, r2);
            Sperp2 = exp(-Lperp2.*Gperp2_mag2);           
            
            % compute signal as the product of the value along and
            % perpendicular to the cylinder axis
            gf_gs0 = Spar.*Sperp1.*Sperp2;
            sig = s0*gf_gs0;
            
            if nargout==2
                out.gf_gs0 = gf_gs0;
                out.U = U;
                out.G = G;
                out.Gpar_mag2 = Gpar_mag2;
                out.Gperp1_mag2 = Gperp1_mag2;
                out.Gperp2_mag2 = Gperp2_mag2;
                out.Lpar = Lpar;
                out.Lperp1 = Lperp1;
                out.Lperp2 = Lperp2;
                out.nom1 = nom1;
                out.denom1 = denom1;
                out.nom2 = nom2;
                out.denom2 = denom2;
            end
        end % of synthesize
        
        function jac = jacobian(obj, params, schemefile)
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
            
            [f0, out] = obj.synthesize(params, schemefile);
            
            % read params into individual model parameters
            s0      = params(1);   % b0 signal
            diffPar = params(2);   % diffusivity [s/m^2]
            r1      = params(3);   % radius along the long axis [m]
            r2      = params(4);   % radius along the short axis [m]
            theta   = params(5);   % elevation angle [radian]
            phi     = params(6);   % azimuth angel [radian]
            alpha   = params(7);   % plane angle [radian]
            
            % read intermediate variables
            gf_gs0 = out.gf_gs0;
            U = out.U;
            G = out.G;
            Gpar_mag2 = out.Gpar_mag2;
            Gperp1_mag2 = out.Gperp1_mag2;
            Gperp2_mag2 = out.Gperp2_mag2;
            Lpar = out.Lpar;
            Lperp1 = out.Lperp1;
            Lperp2 = out.Lperp2;
            nom1 = out.nom1;
            denom1 = out.denom1;
            nom2 = out.nom2;
            denom2 = out.denom2;
            
            % initialise the jac with zeros
            jac = zeros(schemefile.measurementsNum(), obj.nParams);
            
            % gradient wrt to s0
            jac(:,1) = gf_gs0;

            % gradient wrt diffPar  
            gLpar_gDiff = cylinder.get_gLpar_gDiff(schemefile);
            gLperp1_gDiff = cylinder.get_gLperp_gDiff(schemefile, diffPar, r1, nom1, denom1);
            gLperp2_gDiff = cylinder.get_gLperp_gDiff(schemefile, diffPar, r2, nom2, denom2);    
            gf_gDiff = -(gLpar_gDiff.*Gpar_mag2 + gLperp1_gDiff.*Gperp1_mag2 + gLperp2_gDiff.*Gperp2_mag2).*f0;
            jac(:,2) = gf_gDiff;
            

            % gradient wrt r1
            gLperp1_gR1 = cylinder.get_gLperp_gR(schemefile, diffPar, r1, nom1, denom1);
            gf_gR1 = -(gLperp1_gR1.*Gperp1_mag2).*f0;
            jac(:,3) = gf_gR1;

            % gradient wrt r2
            gLperp2_gR2 = cylinder.get_gLperp_gR(schemefile, diffPar, r2, nom2, denom2);
            gf_gR2 = -(gLperp2_gR2.*Gperp2_mag2).*f0;
            jac(:,4) = gf_gR2;
            
            
            
            % compute gradients of V1, V2, and V3 wrt theta, phi, alpha
            [gU1, gU2, gU3] = math.getOrientationJacobian(theta, phi, alpha);
            gf_gU1 = repmat(-2*Lpar.*(G*U(:,1)).*f0,1,3).*G;
            gf_gU2 = repmat(-2*Lperp1.*(G*U(:,2)).*f0,1,3).*G;
            gf_gU3 = repmat(-2*Lperp2.*(G*U(:,3)).*f0,1,3).*G;
            
            % gradient wrt theta
            jac(:,5) = gf_gU1*gU1(:,1) + gf_gU2*gU2(:,1) + gf_gU3*gU3(:,1);
            % gradient wrt phi
            jac(:,6) = gf_gU1*gU1(:,2) + gf_gU2*gU2(:,2) + gf_gU3*gU3(:,2);
            % gradient wrt alpha
            jac(:,7) = gf_gU1*gU1(:,3) + gf_gU2*gU2(:,3) + gf_gU3*gU3(:,3);
                           
        end % of jacobian
    end% of methods (public)
    %\\
end % of CylinderECS class
