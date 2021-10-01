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
        nHyperparams = 1;
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
            obj.hyperparams = 1202;
            obj.hyperparamsName = {'npts'};
        end
        
        function [s, tmpAccessMemory] = synthesize(obj, params, schemefile)
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
            
            assert(isa(schemefile, 'scheme'), ...
                'MATLAB:tensor:invalidInputArgument',...
                'Scheme file should be of type scheme.');
            
            assert(strcmp(schemefile.type, 'stejskal-tanner'), ...
                'MATLAB:tensor:invalidInputArgument',...
                'Scheme file should be of type stejskal-tanner.');
            
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
            npts = obj.hyperparams(1);
           
            [orientations, weights] = math.sampleSphere('lebedev', npts);            
            probs = bingham.pdf(...
                [theta, phi, alpha, 0, kappa1, kappa2]', orientations, weights);
            
            
            % gradient matrix [nScheme x 3] for each measurment
            G_dir = schemefile.ghat; 
            G_mag = schemefile.gmag;
            
            % comput the gradient along and perpendicular
            % to the cylinder axis [nScheme x nInstance]
            Gpar_mag = (G_dir*orientations').*G_mag;
            Gpar_mag2 = Gpar_mag.^2; 
            Gperp_mag2 = G_mag.^2-Gpar_mag2; 
            
            % compute the signal parallel to the cylinder axis
            % [nScheme x nInstance]
            Lpar = cylinder.get_Lpar(schemefile, diffPar);
            sPar = exp(-Lpar.*Gpar_mag2);
            
            % compute the signal perpendicular to the cylinder axis
            % [1 x length(Jp1ROOTS)]
            [Lperp, nom, denom] = cylinder.get_Lperp(schemefile, diffPar, r);
            sPerp = exp(-Lperp.*Gperp_mag2);
            
            % compute signal as the product of the value along and
            % perpendicular to the cylinder axis
            sig = s0 * sPerp .* sPar;
            
            % compute the average wrt Bingham Distribution
            %weights = probs.*sin(thisTheta(:));
            s = sig*probs;
            gf0_gs0 = (sPerp.*sPar)*probs;
            
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
                tmpAccessMemory.orientations = orientations;
                tmpAccessMemory.weights = weights;
            end
        end
        
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
            %                   jac: Numerical matrix of size M x nParams.
            %
            % synthesize signal
            [~, tmpAccessMemory] = obj.synthesize(params, schemefile);
            
            % read intermediate variables
            gf0_gs0 = tmpAccessMemory.gf0_gs0;
            nom = tmpAccessMemory.nom;
            denom = tmpAccessMemory.denom;
            Gpar_mag2 = tmpAccessMemory.Gpar_mag2;
            Gperp_mag2 = tmpAccessMemory.Gperp_mag2;
            sig = tmpAccessMemory.sig;
            probs = tmpAccessMemory.probs;
            orientations = tmpAccessMemory.orientations;
            weights = tmpAccessMemory.weights;
            
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
            jac = zeros(schemefile.measurementsNum(), obj.nParams);
            
            %%
            % gSgS0 : gradient wrt s0
            jac(:,1) = gf0_gs0;
            
            % gradient wrt diffPar
            % first compute gradients for each pair f(theta', phi') and then
            % average over the Bingham distriubtion.
            gLpar_gDiffPar = cylinder.get_gLpar_gDiff(schemefile);
            gLperp_gDiffPar = cylinder.get_gLperp_gDiff(schemefile, diffPar, r, nom, denom);
            
            % gradient wrt diffPar      
            gf_gDiffPar = -(gLpar_gDiffPar.*Gpar_mag2 + gLperp_gDiffPar.*Gperp_mag2).*sig;
            
            % gSgDiffPar: average gradient over the Bingham distribution
            jac(:,2) = gf_gDiffPar*probs; 
            
            % gradient wrt r
            % first compute gradients for each pair f(theta', phi') and then
            % average over the Bingham distriubtion.
            gLperp_gR = cylinder.get_gLperp_gR(schemefile, diffPar, r, nom, denom);
            gf_gR = -(gLperp_gR.*Gperp_mag2).*sig;
            
            % gSgR
            jac(:,3)  = gf_gR*probs;
            
            %% gradient wrt Bingham distribution parameters
            gBingham = bingham.jacobian([theta, phi, alpha, 0, kappa1, kappa2]',...
                orientations, weights);
               
            % gSgTheta 
            jac(:,4)   = sig*gBingham(:,1); 
            
            % gSgPhi
            jac(:,5)   = sig*gBingham(:,2);
                     
            % gSgAlpha
            jac(:,6)   = sig*gBingham(:,3); 
            
            % gSgKappa1
            jac(:,7)   = sig*gBingham(:,5);
            
            % gSgKappa2
            jac(:,8)   = sig*gBingham(:,6);
        end
        
    end%of methods (public)
    %\\
end%of cylinderBDA class

