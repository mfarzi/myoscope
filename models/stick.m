classdef stick < compartment
    % STICK 
    % Cylinder Model with Zero Radius
    %
    %   A STICK object is a basic COMPARTMENT object represenitng particles
    %   diffusion in a cylinder with zero radius. This class implements
    %   inherited abstract methods and properties. This class should be 
    %   used as input to MULTICOMPARTMENT class.
    %
    %   Model Parameters:
    %       name              - class name
    %       s0                - b0-signal with no diffusion weight or
    %                           signal attenuation at b-val=0 
    %       diff              - diffusivity along the cylinder axis [m^2/s]
    %       theta             - elevation angle; the angle between the
    %                           cylinder axis (n1) and the z-axis. [radian]
    %       phi               - azimuth angle; the angle between the x-axis
    %                           and the n1 projection onto the xy-plane. 
    %                           [radian]
    %       fitter            - an "optimizer" object to fit model
    %                           parameters
    %
    %   For mathematical background see
    %       Panagiotaki, E., Schneider, T., Siow, B., Hall, M.G., Lythgoe,
    %       M.F. and Alexander, D.C., "Compartment models of the diffusion
    %       MR signal in brain white matter: a taxonomy and comparison.",
    %       Neuroimage, 59(3), pp.2241-2254, 2012.
    %
    %   See also: multicompartment, cylinder, cylinderECS, cylinderGDR,
    %             cylinderBDA
    %
    % Mohsen Farzi
    % Email: m.farzi@leeds.ac.uk
    
    properties (SetAccess='protected')
        name = 'stick';    
        nCompartments = 1;     % number of basic COMPARTMETNT objects
        nParams = 4;
        nHyperparams = 0;
    end
     
    properties (Access = 'protected')
        comp = [];             % List of COMPARTMENT object 
        links;
        hyperparams = [];
        hyperparamsName = {};
    end     
    
    methods
        function obj = stick()
            %STICK Construct Function.
            %
            %   stick() constructs a single-compartment model
            
            % initialise the linking functions to map constrained model
            % parameters to unconstrained optimisation variables.
            s0    = linker('s0'  , 'bounded', 0    , 1);
            diff  = linker('diff', 'bounded', 1e-10, 3e-9);
            theta = linker('theta');
            phi   = linker('phi');
            obj.links = [s0; diff; theta; phi];
        end
        
        function s = synthesize(obj, params, scheme)
            % synthesize(params, scheme, hyperparams) return DW-MR signals
            % from a cylinder with r=0. Diffusion is only allowed 
            % along the cylinder axis; no signal attenuation occurs in 
            % plane perpendicular to the cylinder axis.
            %       Input arguments:
            %                params: Numerical column vector of all model
            %                        parameters
            %                scheme: Diffusion scheme
            %
            %      output arguments:
            %                     s: Numerical column vector of synthesied
            %                        signal
            
            % validate inputs
            validateattributes(params, {'double'},...
                {'column', 'nrows', obj.nParams},...
                'stick.synthesize', 'params');
            
            % read params into individual model parameters
            s0    = params(1);   % normalised b0 signal
            diff  = params(2);   % diffusivity [s/m^2]
            theta = params(3);   % elevation angle [radian]
            phi   = params(4);   % azimuth angel [radian]
           
            % compute the cylinderal axis
            N = math.getUnitFrame(theta, phi, 0);
            n = N(:,1);

            % compute the signal parallel to the cylinder axis
            bPar = (scheme.DELTA-scheme.delta/3).*(math.GAMMA*scheme.delta).^2*diff;
            G = [scheme.x, scheme.y, scheme.z].*scheme.G_mag;
            s = s0*exp(-bPar.*(G*n).^2);
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
            
            % validate inputs
            validateattributes(params, {'double'},...
                {'column', 'nrows', obj.nParams},...
                'stick.synthesize', 'params');
            
            % read params into individual model parameters
            s0    = params(1);   % normalised b0 signal
            diff  = params(2);   % diffusivity [s/m^2]
            theta = params(3);   % elevation angle [radian]
            phi   = params(4);   % azimuth angel [radian]
            
            % initialise the jac with zeros
            jac = zeros(size(scheme, 1), obj.nParams);
            
            % compute the cylinderal axis
            N = math.getUnitFrame(theta, phi, 0);
            n = N(:,1);
            
            % compute the signal parallel to the cylinder axis
            bPar = (scheme.DELTA-scheme.delta/3).*(math.GAMMA*scheme.delta).^2;
            G = [scheme.x, scheme.y, scheme.z].*scheme.G_mag;
            f0 = s0*exp(-diff*bPar.*(G*n).^2);
            
            
            % gradient wrt s0
            jac(:,1) = exp(-diff*bPar.*(G*n).^2);
            
            % gradient wrt diff 
            gf_gDiff = (-bPar.*(G*n).^2).*f0;
            jac(:,2) = gf_gDiff;
            
            % gradient wrt theta
            gf_gN = repmat(-2*diff*bPar.*(G*n).*f0,1,3).*G;
            
            gN_gTheta = [cos(phi)*cos(theta); ...
                         sin(phi)*cos(theta); ...
                        -sin(theta)];
            jac(:,3) = gf_gN*gN_gTheta;
            
            % gradient wrt phi
            gN_gPhi   = [-sin(phi)*sin(theta); ...
                          cos(phi)*sin(theta); ...
                           0];
            jac(:,4) = gf_gN * gN_gPhi;           
        end % of jacobian

    end%of methods (public)
end
