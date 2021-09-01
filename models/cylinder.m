classdef cylinder < compartment
    % CYLINDER 
    % (Cylinder Model with Gaussian Phase Distribution)
    %
    %   A CYLINDER object is a basic COMPARTMENT object represenitng 
    %   particles diffusion in a cylinder. This class implements inherited
    %   abstract methods and properties. This class should be used as input
    %   to MULTICOMPARTMENT class.
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
    %
    %   For mathematical background see
    %       Vangelderen, P., DesPres, D., Vanzijl, P.C.M. and Moonen, C.,
    %       "Evaluation of restricted diffusion in cylinders.
    %       Phosphocreatine in rabbit leg muscle.", Journal of Magnetic
    %       Resonance, Series B, 103(3), pp.255-260, 1994.
    %
    %   See also: multicompartment, cylinderECS, cylinderBDA, cylinderGDR,
    %             stick
    %
    % Mohsen Farzi
    % Email: m.farzi@leeds.ac.uk
    
    properties (SetAccess='protected')
        name = 'cylinder';    
        nCompartments = 1;     % number of basic COMPARTMETNT objects 
        nParams = 5;
        nHyperparams = 0;
    end
    
    properties (Access = 'protected')
        comp = [];             % List of COMPARTMENT object 
        links;
        hyperparams = [];
        hyperparamsName = {};
    end
    
    methods
        function obj = cylinder()
            %CYLINDER Construct Function.
            %
            %   cylinder() constructs a single-compartment model
            
            % initialise the linking functions to map constrained model
            % parameters to unconstrained optimisation variables.
            s0        = linker('s0'       , 'bounded', 0    , 1);
            diffPar   = linker('diffPar'  , 'bounded', 1e-10, 3e-9);
            r         = linker('r'        , 'bounded', 1e-6 , 20e-6);
            theta     = linker('theta');
            phi       = linker('phi');            
            obj.links = [s0; diffPar; r; theta; phi];
        end
        
        function [sig, out] = synthesize(obj, params, scheme)
            % synthesize(params, scheme, hyperparams) return DW-MR signals
            % from a cylinder with Gaussian Phase Distribution. The signal 
            % is modeled as the product of the signals parallel and 
            % perpendicular to the cylinder axis.
            %       Input arguments:
            %                params: Numerical column vector of all model
            %                        parameters
            %                scheme: Diffusion scheme
            %
            %      output arguments:
            %                     s: Numerical column vector of synthesied
            %                        signal
            %                   out: store intermediate variables to be
            %                        used with method 'jacobian'
            %  
            
            % validate inputs
            validateattributes(params, {'double'},...
                {'column', 'nrows', obj.nParams},...
                'cylinder.synthesize', 'params');
            
            % read params into individual model parameters
            s0      = params(1);   % b0 signal
            diffPar = params(2);   % diffusivity [s/m^2]
            r       = params(3);   % radius [m]
            theta   = params(4);   % elevation angle
            phi     = params(5);   % azimuth angel
            
            % gradient matrix [nMeasurement x 3] 
            G = [scheme.x, scheme.y, scheme.z].*scheme.G_mag;
            
            % compute the cylinderal axis
            n = math.get_u1(theta, phi);
            
            % comput the gradient along and perpendicular
            % to the cylinder axis 
            Gpar_mag2 = (G*n).^2; 
            Gperp_mag2 = sum(G.^2,2)-Gpar_mag2;
            
            % compute the signal parallel to the cylinder axis
            Lpar = cylinder.get_Lpar(scheme, diffPar);
            Spar = exp(-Lpar.*Gpar_mag2);

            % compute the signal perpendicular to the cylinder axis
            [Lperp, nom, denom] = cylinder.get_Lperp(scheme, diffPar, r);
            Sperp = exp(-Lperp.*Gperp_mag2);
            
            % compute signal as the product of the value along and
            % perpendicular to the cylinder axis
            s = Sperp.*Spar;
            sig = s0 * s;
            
            if nargout==2
                out.s = s;
                out.G = G;
                out.Gpar_mag2 = Gpar_mag2;
                out.Gperp_mag2 = Gperp_mag2;
                out.Lpar = Lpar;
                out.Lperp = Lperp;
                out.nom = nom;
                out.denom = denom;
                out.n = n;
            end
        end % of synthesize
        
        function [jac, f0] = jacobian(obj, params, scheme)
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
            
            % run synthesize
            [f0, out] = obj.synthesize(params, scheme, []);
            
            % read intermediate variables
            gf0_gs0 = out.s;
            G = out.G;
            Gpar_mag2 = out.Gpar_mag2;
            Gperp_mag2 = out.Gperp_mag2;
            Lpar = out.Lpar;
            Lperp = out.Lperp;
            nom = out.nom;
            denom = out.denom;
            n = out.n;
            
            % read params into individual model parameters
            s0      = params(1);   % b0 signal
            diffPar = params(2);   % diffusivity [s/m^2]
            r       = params(3);   % radius [m]
            theta   = params(4);   % elevation angle
            phi     = params(5);   % azimuth angel
            
            % initialise the jac with zeros
            jac = zeros(size(scheme, 1), obj.nParams);            
            
            % gradient wrt to s0
            jac(:,1) = gf0_gs0;
            
            % gLpar_gDiff
            gLpar_gDiff = cylinder.get_gLpar_gDiff(scheme);
     
            % gLperp_gDiff
            gLperp_gDiff = cylinder.get_gLperp_gDiff(scheme, diffPar, r, nom, denom);
            
            % gradient wrt diffPar      
            gf_gDiff = -(gLpar_gDiff.*Gpar_mag2 + gLperp_gDiff.*Gperp_mag2).*f0;
            jac(:,2) = gf_gDiff;
            
            % gLperp_gR
            gLperp_gR = cylinder.get_gLperp_gR(scheme, diffPar, r, nom, denom);
            
            % gradient wrt R
            gf_gR = -(gLperp_gR.*Gperp_mag2).*f0;
            jac(:,3) = gf_gR;
            
            % gradient of f wrt N
            gf_gN = (repmat(2*(Lperp-Lpar).*(G*n).*f0,1,3).*G);
            
            % gradient wrt to theta
            gN_gTheta = math.get_gU1_gTheta(theta, phi);
            jac(:,4) = gf_gN*gN_gTheta;         
            
            % Gradient wrt to phi
            gN_gPhi = math.get_gU1_gPhi(theta, phi);
            jac(:,5) = gf_gN*gN_gPhi;                   
        end % of jacobian
    end%of method (public)
    
    methods (Static)
        function Lpar = get_Lpar(scheme, diffPar)
            %
            Lpar = (scheme.DELTA-scheme.delta/3).*(math.GAMMA*scheme.delta).^2*diffPar;
        end
        
        function g = get_gLpar_gDiff(scheme)
            %
            g = (scheme.DELTA-scheme.delta/3).*(math.GAMMA*scheme.delta).^2;
        end
        
        function [Lperp, nom, denom] = get_Lperp(scheme, diffPar, r)
            %
            beta = math.Jp1ROOTS'/r;          
            B2   = beta.^2;
            
            %[nScheme x length(Jp1ROOTS)]
            diff_dot_delta_dot_beta2 = diffPar*(scheme.delta*B2);
            diff_dot_DELTA_dot_beta2 = diffPar*(scheme.DELTA*B2);
            
            nom = 2*diff_dot_delta_dot_beta2 - 2  ...
                + 2*exp(-diff_dot_delta_dot_beta2) ...
                + 2*exp(-diff_dot_DELTA_dot_beta2) ...
                - exp(-diff_dot_DELTA_dot_beta2+diff_dot_delta_dot_beta2) ...
                - exp(-diff_dot_DELTA_dot_beta2-diff_dot_delta_dot_beta2);
            
            % [length(Jp1ROOTS) x 1]
            denom = diffPar^2 * ((beta.^6).*(math.Jp1ROOTS.^2-1)');
            
            %[nScheme x 1]
            Lperp = 2*math.GAMMA^2*sum(nom./denom, 2);
            %Lperp = 2*math.GAMMA^2*nom*(1./denom)';
        end
        
        function gLperp_gDiff = get_gLperp_gDiff(scheme, diffPar, r, nom, denom)
            %
            beta = math.Jp1ROOTS'/r;          
            B2   = beta.^2;
            
            %[nScheme x length(Jp1ROOTS)]
            delta_dot_beta2 = scheme.delta*B2;
            DELTA_dot_beta2 = scheme.DELTA*B2;
            
            gNom_gDiffPar = 2*delta_dot_beta2 ...
                          - 2*exp(-diffPar*delta_dot_beta2).*delta_dot_beta2 ...
                          - 2*exp(-diffPar*DELTA_dot_beta2).*DELTA_dot_beta2 ...
                          +   exp(-diffPar*(DELTA_dot_beta2-delta_dot_beta2)).*(DELTA_dot_beta2-delta_dot_beta2)...
                          +   exp(-diffPar*(DELTA_dot_beta2+delta_dot_beta2)).*(DELTA_dot_beta2+delta_dot_beta2);
            
            gDenom_gDiffPar = 2*diffPar*beta.^6.*(math.Jp1ROOTS.^2-1)';
            tmp = (gNom_gDiffPar.*denom - nom.*gDenom_gDiffPar)./(denom.^2);
            gLperp_gDiff = 2*math.GAMMA^2*sum(tmp, 2);
        end
        
        function gLperp_gR = get_gLperp_gR(scheme, diffPar, r, nom, denom)
            beta = math.Jp1ROOTS'/r;          
            B2   = beta.^2;
            
            %[nScheme x length(Jp1ROOTS)]
            delta_dot_beta2 = scheme.delta*B2;
            delta_dot_beta  = scheme.delta*beta;
            DELTA_dot_beta2 = scheme.DELTA*B2;
            DELTA_dot_beta  = scheme.DELTA*beta;
            
            gNom_gBm = 4*diffPar*delta_dot_beta ...
                     - 2*exp(-diffPar*delta_dot_beta2).*(2*diffPar*delta_dot_beta) ...
                     - 2*exp(-diffPar*DELTA_dot_beta2).*(2*diffPar*DELTA_dot_beta) ...
                     +   exp(-diffPar*(DELTA_dot_beta2-delta_dot_beta2)).*(2*diffPar*(DELTA_dot_beta-delta_dot_beta)) ...
                     +   exp(-diffPar*(DELTA_dot_beta2+delta_dot_beta2)).*(2*diffPar*(DELTA_dot_beta+delta_dot_beta));
            
            gDenom_gBm = 6*diffPar^2*(beta.^5.*(math.Jp1ROOTS.^2-1)');
            
            gLperp_gBm = 2 * math.GAMMA^2 * ...
                        (gNom_gBm.*denom - nom .* gDenom_gBm)./(denom.^2);
            
            % gLperp_gR
            gBm_gR = -math.Jp1ROOTS/r^2;
            gLperp_gR = gLperp_gBm*gBm_gR;
        end
        
    end % of methods (protected)
    %\\
end%of class

