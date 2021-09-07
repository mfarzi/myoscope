classdef tensor < compartment
    % TENSOR 
    %
    %   A TENSOR object is a basic COMPARTMENT object represenitng 
    %   particles diffusion in a medium using a diffusion tensor (DT). This
    %   class implements inherited abstract methods and properties. This
    %   class should be used as input to MULTICOMPARTMENT class.
    %
    %   Model Parameters:
    %       s0                - Normalised b0-signal 
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
    %
    %   For mathematical background see
    %       Basser, P.J., Mattiello, J. and LeBihan, D., "MR diffusion 
    %       tensor spectroscopy and imaging.", Biophysical journal, 66(1),
    %       pp.259-267, 1994.
    %
    %
    %   See also: compartment, ball, zeppelin, pancake, multicompartment
    %
    % Mohsen Farzi
    % Email: m.farzi@leeds.ac.uk   
    
    properties (SetAccess = 'protected')
        name = 'tensor';  
        nCompartments = 1;     % number of basic COMPARTMETNT objects 
        nParams = 7;           % Number of model parameters
        nHyperparams = 0;      % Number of model hyper-parameters
    end
    
    properties (Access = 'protected')
        comp = [];             % List of COMPARTMENT object 
        links;
        hyperparams = [];
        hyperparamsName = {};
    end
    
    methods
        function obj = tensor()
            %TENSOR Construct Function.
            %
            %   tensor() constructs a single-compartment model
            
            % initialise the linking functions to map constrained model
            % parameters to unconstrained optimisation variables.
            s0        = linker('s0'       , 'bounded', 0    , 1);
            diffPar   = linker('diffPar'  , 'bounded', 1e-10, 3e-9);
            diffPerp1 = linker('diffPerp1', 'bounded', 1e-10, 3e-9);
            diffPerp2 = linker('diffPerp2', 'bounded', 1e-10, 3e-9);
            theta     = linker('theta');
            phi       = linker('phi');
            alpha     = linker('alpha');
            links = [s0; diffPar; diffPerp1; diffPerp2; theta; phi; alpha];
            links.addConstraint('diffPar>=diffPerp1');
            links.addConstraint('diffPerp1>=diffPerp2');
            obj.links = links;
        end

        function [sig, out] = synthesize(obj, params, scheme)
            % synthesize(params, scheme, hyperparams) return DW-MR signals
            % from a single diffusion tensor.
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
                'tensor.synthesize', 'params');
            
            % read params into individual model parameters
            s0          = params(1);   % b0 signal
            diffPar     = params(2);   % diffusivity [s/m^2]
            diffPerp1   = params(3);   % diffusivity [s/m^2]
            diffPerp2   = params(4);   % diffusivity [s/m^2]
            theta       = params(5);   % elevation angle
            phi         = params(6);   % azimuth angel
            alpha       = params(7);   % angele beween v1_zrot and v2
            
            G_dir = [scheme.x, scheme.y, scheme.z];
            
            b = (scheme.DELTA-scheme.delta/3).* ...
                ((scheme.delta .*scheme.G_mag)*math.GAMMA).^2;
            
            % compute the orhtonormal basis
            U = math.getUnitFrame(theta, phi, alpha);
            
            gf_gs0 = exp(-b.*(diffPar  *(G_dir*U(:,1)).^2 + ...
                         diffPerp1*(G_dir*U(:,2)).^2 + ...
                         diffPerp2*(G_dir*U(:,3)).^2)); 
            
            sig = s0*gf_gs0;
            
            if nargout==2
                out.gf_gs0 = gf_gs0;
                out.U = U;
                out.G_dir = G_dir;
                out.b = b;
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
            
            [f0, out] = obj.synthesize(params, scheme);
            
            % read params
            s0          = params(1);   % b0 signal
            diffPar     = params(2);   % diffusivity [s/m^2]
            diffPerp1   = params(3);   % diffusivity [s/m^2]
            diffPerp2   = params(4);   % diffusivity [s/m^2]
            theta       = params(5);   % elevation angle
            phi         = params(6);   % azimuth angel
            alpha       = params(7);   % angele beween v1_zrot and v2
            
            % read out variables
            gf_gs0 = out.gf_gs0;
            U = out.U;
            G_dir = out.G_dir;
            b = out.b;
            
            % initialise the jac with zeros
            jac = zeros(size(scheme,1), obj.nParams);
                      
            % gradient wrt to s0
            jac(:,1) = gf_gs0;
            
            % gradient wrt diffPar
            gf_gDiffPar = (-b.*((G_dir*U(:,1)).^2)).*f0;
            jac(:,2) = gf_gDiffPar;
            
            % gradient wrt diffPerp1
            gf_gDiffPerp1 = (-b.*((G_dir*U(:,2)).^2)).*f0;
            jac(:,3) = gf_gDiffPerp1;
            
            % gradient wrt diffPerp2
            gf_gDiffPerp2 = (-b.*((G_dir*U(:,3)).^2)).*f0;
            jac(:,4) = gf_gDiffPerp2;
            
            % gradient wrt theta, phi, and alpha
            % using chain rule, comput gf_gV1, gf_gV2, gf_gV3.
            gf_gV1 = repmat(-2*diffPar*b.*(G_dir*U(:,1)).*f0,1,3).*G_dir;
            gf_gV2 = repmat(-2*diffPerp1*b.*(G_dir*U(:,2)).*f0,1,3).*G_dir;
            gf_gV3 = repmat(-2*diffPerp2*b.*(G_dir*U(:,3)).*f0,1,3).*G_dir;
            
            % compute gradients of V1, V2, and V3 wrt theta, phi, alpha
            [gV1, gV2, gV3] = math.getOrientationJacobian(theta, phi, alpha);
            
            % gradient wrt theta           
            jac(:,5) = gf_gV1*gV1(:,1) + gf_gV2*gV2(:,1) + gf_gV3*gV3(:,1);          
            
            % gradient wrt phi         
            jac(:,6) = gf_gV1*gV1(:,2) + gf_gV2*gV2(:,2) + gf_gV3*gV3(:,2);
            
            % gradient wrt alpha
            jac(:,7) = gf_gV1*gV1(:,3) + gf_gV2*gV2(:,3) + gf_gV3*gV3(:,3);
           
        end % of jacobian
    end % of methods (public)
end
