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
            diffPar  = linker('diffPar', 'bounded', 1e-10, 3e-9);
            theta = linker('theta');
            phi   = linker('phi');
            obj.links = [s0; diffPar; theta; phi];
        end
        
        function [sig, tmpAccessMemory] = synthesize(obj, params, schemefile)
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
            
            assert(isa(schemefile, 'scheme'), ...
                'MATLAB:stick:invalidInputArgument',...
                'Scheme file should be of type scheme.');
            
            % read params into individual model parameters
            s0       = params(1);   % normalised b0 signal
            diffPar  = params(2);   % diffusivity [s/m^2]
            theta    = params(3);   % elevation angle [radian]
            phi      = params(4);   % azimuth angel [radian]
            
            % read scheme file
            b = schemefile.bval;
            ghat = schemefile.ghat;
            
            % compute the unit frame
            U = math.getUnitFrame(theta, phi, 0);

            % compute the signal parallel to the primary axis
            gf_gs0 = exp(-b.*(diffPar  *(ghat*U(:,1)).^2));
            sig = s0*gf_gs0;
            if nargout==2
                tmpAccessMemory.gf_gs0 = gf_gs0;
                tmpAccessMemory.U = U;
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
            
            [f0, out] = obj.synthesize(params, schemefile);
            
            % read params into individual model parameters
            s0      = params(1);   % normalised b0 signal
            diffPar = params(2);   % diffusivity [s/m^2]
            theta   = params(3);   % elevation angle [radian]
            phi     = params(4);   % azimuth angel [radian]
            
            % read scheme file
            ghat = schemefile.ghat;
            b = schemefile.bval;
            
            % read out variables
            gf_gs0 = out.gf_gs0;
            U = out.U;
            
            % initialise the jac with zeros
            jac = zeros(schemefile.measurementsNum(), obj.nParams);
            
            % gradient wrt s0
            jac(:,1) = gf_gs0;
            
            % gradient wrt diffPar 
            gf_gDiffPar = (-b.*((ghat*U(:,1)).^2)).*f0;
            jac(:,2) = gf_gDiffPar;
            
            % gradient wrt theta and phi
            gf_gU1 = repmat(-2*diffPar*b.*(ghat*U(:,1)).*f0,1,3).*ghat;
            gU1    = math.getOrientationJacobian(theta, phi, 0);
            
            jac(:,3) = gf_gU1*gU1(:,1);
            jac(:,4) = gf_gU1*gU1(:,2);           
        end % of jacobian

    end%of methods (public)
end
