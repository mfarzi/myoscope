classdef ball < compartment
    % BALL
    %
    %   A BALL object is a basic COMPARTMENT object represenitng particles 
    %   diffusion in an isotropic, unconstrained medium. This class
    %   implements inherited abstract methods and properties. This class
    %   should be used as an input to the MULTICOMPARTMENT class.
    %
    %   For mathematical background see
    %       Panagiotaki, E., Schneider, T., Siow, B., Hall, M.G., Lythgoe,
    %       M.F. and Alexander, D.C., "Compartment models of the diffusion
    %       MR signal in brain white matter: a taxonomy and comparison.",
    %       Neuroimage, 59(3), pp.2241-2254, 2012.
    %
    %   Model Parameters:
    %       s0        - Normalised b0-signal 
    %       diff      - diffusivity of ball [m^2/s]
    %
    %   See also: multicompartment, tensor, zeppelin, pancake
    %
    % Mohsen Farzi
    % Email: m.farzi@leeds.ac.uk
    
    properties (SetAccess='protected')
        name = 'ball';    
        nCompartments = 1;     % number of basic COMPARTMETNT objects 
        nParams = 2;
        nHyperparams = 0;
    end
    
    properties (Access = 'protected')
        comp = [];             % List of COMPARTMENT object 
        links;
        hyperparams = [];
        hyperparamsName = {};
    end
    
    methods
        function obj = ball()
            %BALL Construct Function.
            %
            %   ball() constructs a single-compartment model
            
            % initialise the linking functions to map constrained model
            % parameters to unconstrained optimisation variables.
            s0    = linker('s0'  , 'bounded', 0    , 1);
            diff  = linker('diff', 'bounded', 1e-10, 3e-9);
            links = [s0; diff];
            obj.links = links;
        end
        
        function s = synthesize(obj, params, schemefile)
            % synthesize(params, scheme, hyperparams) return DW-MR signals
            % from a single isotropic diffusion tensor.
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
                'ball.synthesize', 'params');
            
            assert(isa(schemefile, 'scheme'), ...
                'MATLAB:ball:invalidInputArgument',...
                'Scheme file should be of type scheme.');
            
            % read params into individual model parameters
            s0   = params(1);   % b0 signal
            diff = params(2);   % diffusivity [s/m^2]
            
            b = schemefile.bval;
            s = s0*exp(-b*diff);
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
            %                   jac: Numerical matrix of size M x nParams.

            % validate inputs
            validateattributes(params, {'double'},...
                {'column', 'nrows', obj.nParams},...
                'ball.synthesize', 'params');
            
            assert(isa(schemefile, 'scheme'), ...
                'MATLAB:tensor:invalidInputArgument',...
                'Scheme file should be of type scheme.');
            
            % read params into individual model parameters
            s0   = params(1);   % b0 signal
            diff = params(2);   % diffusivity [s/m^2]
            
            % initialise the jac with zeros
            jac = zeros(schemefile.measurementsNum(), obj.nParams);
            
            b = schemefile.bval;
            
            % gradient wrt s0
            gf_gs0 = exp(-b*diff);
            jac(:,1) = gf_gs0;
            
            % gradient wrt diff  
            f0 = s0*gf_gs0;
            gf_gDiff = -b.*f0;
            jac(:,2) = gf_gDiff;
            %\\
        end%of jacobian
        
    end%of methods (public) 
end%of BALL class
