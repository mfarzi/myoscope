classdef stickBDA < compartment
    % TENSOR 
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
        name = 'stickBDA';    
        nCompartments = 1;     % number of basic COMPARTMETNT objects 
        nParams = 7;
        nHyperparams = 1;
    end
    
    properties (Access = 'protected')
        comp = [];             % List of COMPARTMENT object 
        links;
        hyperparams = [];
        hyperparamsName = {};
    end
    
    methods 
        function obj = stickBDA(varargin) 
            %STICKBDA Construct Function.
            %
            %   stickBDA() constructs a single-compartment model
            
            % initialise the linking functions to map constrained model
            % parameters to unconstrained optimisation variables.
            s0        = linker('s0'       , 'bounded', 0    , 1);
            diffPar   = linker('diffPar'  , 'bounded', 1e-10, 3e-9);
            theta     = linker('theta');
            phi       = linker('phi');
            alpha     = linker('alpha');
            kappa1  = linker('kappa1' , 'bounded', 0, 128); 
            kappa2  = linker('kappa2' , 'bounded', 0, 128);
            obj.links = [s0; diffPar; theta; phi; alpha; kappa1; kappa2];
            obj.links.addConstraint('kappa2>=kappa1');             
            
            % set hparams and hparams names
            obj.hyperparams = 1202;
            obj.hyperparamsName = {'npts'};
        end%of constructor    
        
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
                'stickBDA.synthesize', 'params');

            % read params into individual model parameters
            s0          = params(1);   % b0 signal
            diffPar     = params(2);   % diffusivity [s/m^2]
            theta       = params(3);   % elevation angle
            phi         = params(4);   % azimuth angel
            alpha       = params(5);   % angele beween v1_zrot and v2
            kappa1      = params(6);
            kappa2      = params(7);
            
            % read hyper-parameters
            npts = obj.hyperparams(1);
           
            [orientations, weights] = math.sampleSphere('lebedev', npts);            
            binghamParams = [theta, phi, alpha, 0, kappa1, kappa2]';
            probs = bingham.pdf(binghamParams, orientations, weights);
                      
            % compute signal along each orientation
            bval = math.GAMMA^2*(scheme.DELTA-scheme.delta/3).*(scheme.delta.*scheme.G_mag).^2;
            G_dir = [scheme.x, scheme.y, scheme.z];
            sig = exp(-diffPar*bval.*(G_dir*orientations').^2);
            
   
            gf_gs0 = sig*probs;
            s = s0*gf_gs0;
            if nargout == 2
                tmpAccessMemory.gf_gs0 = gf_gs0;
                tmpAccessMemory.sig = sig;
                tmpAccessMemory.probs = probs;
                tmpAccessMemory.G_dir = G_dir;
                tmpAccessMemory.bval = bval;
                tmpAccessMemory.orientations = orientations;
                tmpAccessMemory.weights = weights;
            end
        end%of synthesize
        
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
            
            [~, tmpAccessMemory] = obj.synthesize(params, scheme);
           
             % read params into individual model parameters
            s0          = params(1);   % b0 signal
            diffPar     = params(2);   % diffusivity [s/m^2]
            theta       = params(3);   % elevation angle
            phi         = params(4);   % azimuth angel
            alpha       = params(5);   % angele beween v1_zrot and v2
            kappa1      = params(6);
            kappa2      = params(7);
            
            
            % read intermediate variables
            gf_gs0 = tmpAccessMemory.gf_gs0;
            sig = tmpAccessMemory.sig;
            G_dir = tmpAccessMemory.G_dir;
            bval = tmpAccessMemory.bval;
            probs = tmpAccessMemory.probs;
            orientations = tmpAccessMemory.orientations;
            weights = tmpAccessMemory.weights;
            
            % initialise the jac with zeros
            jac = zeros(size(scheme,1), obj.nParams);
                      
            % gradient wrt to s0
            jac(:,1) = gf_gs0;
            
            % gradient wrt diffPar
            gsig_gdiffPar = (-bval.*(G_dir*orientations').^2).*sig;
            jac(:,2) = s0*gsig_gdiffPar*probs;
            
            % gradient wrt Bingham distribution parameters
            gBingham = bingham.jacobian([theta, phi, alpha, 0, kappa1, kappa2]',...
                orientations, weights);
            
            % gradient wrt theta
            jac(:,3) = s0*sig*gBingham(:,1);
            
            % gradient wrt phi
            jac(:,4) = s0*sig*gBingham(:,2);
            
            % gradient wrt alpah
            jac(:,5) = s0*sig*gBingham(:,3);
            
            % gradient wrt kappa1
            jac(:,6) = s0*sig*gBingham(:,5);
            
            % gradient wrt kappa2
            jac(:,7) = s0*sig*gBingham(:,6);     
         end%of jacobian
                
    end % of method (public)
end%of stickBDA class
