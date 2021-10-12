classdef zeppelinCBD < compartment
    %ZEPPELINCBD: zeppelin tensor Coupled with Bingham Distribution
    %
    %
    %   A ZEPPELINCBD object is a basic COMPARTMENT object represenitng 
    %   particles diffusion in a medium using a cylindrically symmetric 
    %   tensor weighted with a Bingham distribution. This class implements
    %   inherited abstract methods and properties.
    %
    %   Model Parameters:
    %       s0        - Normalised b0-signal 
    %       diffPar   - Parallel diffusivity along the primary axis [m^2/s]
    %       diffPerp  - Perpendicular diffusivity [m^2/s]
    %       theta     - Elevation angle; the angle between the first 
    %                   eigenvector (v1) and the z-axis. [radian]
    %       phi       - Azimuth angle; the angle between the x-axis and
    %                   the v1 projection onto the xy-plane.[radian]
    %       alpha     - The angle between the second eigenvector and v1
    %                   rotated by pi/2 around the z-axis. [radian]
    %       kappa1    - Concentration paramter about v2 for the Bingham  
    %                   distribution
    %       kappa2    - Concentration paramter about v3 for the Bingham  
    %                   distribution
    %
    %   For mathematical background see
    %       M. Tariq, T. Schneider, D.C. Alexander, C.A. Gandini
    %       Wheeler-Kingshott, and H. Zhang, "Bingham-NODDI: Mapping
    %       anisotropic orientation dispersion of neurites using diffusion
    %       MRI," in Neuroimage, vol. 133, No. 1, pp. 207?223, 2016.
    %
    %   See also: compartment, zeppelin, zeppelinBDA
    %
    % Mohsen Farzi
    % Email: m.farzi@leeds.ac.uk
    
    properties (SetAccess='protected')
        name = 'zeppelinCBD';    
        nCompartments = 1;     % number of basic COMPARTMETNT objects 
        nParams = 8;
        nHyperparams = 1;
    end
    
    properties (Access = 'protected')
        comp = [];             % List of COMPARTMENT object 
        links;
        hyperparams = [];
        hyperparamsName = {};
    end
    
    methods 
        function obj = zeppelinCBD() 
            %ZEPPELINBDA Construct Function.
            %
            %   zeppelinBDA() constructs a single-compartment model
            
            % initialise the linking functions to map constrained model
            % parameters to unconstrained optimisation variables.
            s0        = linker('s0'       , 'bounded', 0    , 1);
            diffPar   = linker('diffPar'  , 'bounded', 1e-10, 3e-9);
            diffPerp = linker('diffPerp', 'bounded', 1e-10, 3e-9);
            theta     = linker('theta');
            phi       = linker('phi');
            alpha     = linker('alpha');
            kappa1  = linker('kappa1' , 'bounded', 0, 128); 
            kappa2  = linker('kappa2' , 'bounded', 0, 128);
            obj.links = [s0; diffPar; diffPerp;...
                         theta; phi; alpha; kappa1; kappa2];
            obj.links.addConstraint('diffPar>=diffPerp');
            obj.links.addConstraint('kappa2>=kappa1');             
            
            % set hparams and hparams names
            obj.hyperparams = 1202;
            obj.hyperparamsName = {'npts'};
        end%of constructor    
        
        function [sig, tmpAccessMemory] = synthesize(obj, params, schemefile)
            % synthesize(params, scheme, hyperparams) return DW-MR signals.
            %       Input arguments:
            %                params: Numerical column vector of all model
            %                        parameters
            %                scheme: Diffusion scheme
            %
            %      output arguments:
            %                   sig: Numerical column vector of synthesied
            %                        signal
            
            % validate inputs
            validateattributes(params, {'numeric'},...
                {'column', 'nrows', obj.nParams},...
                'cylinderBDA.synthesize', 'params');
            
            assert(isa(schemefile, 'scheme'), ...
                'MATLAB:zeppelinCBD:invalidInputArgument',...
                'Scheme file should be of type scheme.');

            % read params into individual model parameters
            s0          = params(1);   % b0 signal
            diffPar     = params(2);   % diffusivity [s/m^2]
            diffPerp    = params(3);   % diffusivity [s/m^2]
            theta       = params(4);   % elevation angle
            phi         = params(5);   % azimuth angel
            alpha       = params(6);   % angele beween v1_zrot and v2
            kappa1      = params(7);   % concentration parameter
            kappa2      = params(8);   % concentration parameter
            
            cB = bingham.getC(0, kappa1, kappa2);
            gcB_gK1 = bingham.get_gcB_gk1(0, kappa1, kappa2); 
            gcB_gK2 = bingham.get_gcB_gk2(0, kappa1, kappa2);
            
            % Estimate diffusivity for the Coupled Diffusion Tensor
            cdtDiffPar = diffPerp - (diffPar-diffPerp)*(gcB_gK1/cB);
            cdtDiffPerp1 = diffPerp - (diffPar-diffPerp)*(gcB_gK2/cB); 
            cdtDiffPerp2 = diffPar+ 2*diffPerp - cdtDiffPar - cdtDiffPerp1;
            
            cdt = tensor(); % coupled diffusion tensor
            coupuledTensorParams = [s0, cdtDiffPar, cdtDiffPerp1,...
                cdtDiffPerp2, theta, phi, alpha]';
            sig = cdt.synthesize(coupuledTensorParams, schemefile);
            
            if nargout==2
                tmpAccessMemory.cdtDiffPar = cdtDiffPar;
                tmpAccessMemory.cdtDiffPerp1 = cdtDiffPerp1;
                tmpAccessMemory.cdtDiffPerp2 = cdtDiffPerp2;
                tmpAccessMemory.cB = cB;
                tmpAccessMemory.gcB_gK1 = gcB_gK1;
                tmpAccessMemory.gcB_gK2 = gcB_gK2;
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
            %           hyperparams: non-optimisable model parameters
            %
            %      output arguments:
            %                   jac: Numerical matrix of size M x nParams.
            %
            
            [~, tmpAccessMemory] = obj.synthesize(params, schemefile);  
            
            % read params into individual model parameters
            s0          = params(1);   % b0 signal
            diffPar     = params(2);   % diffusivity [s/m^2]
            diffPerp    = params(3);   % diffusivity [s/m^2]
            theta       = params(4);   % elevation angle
            phi         = params(5);   % azimuth angel
            alpha       = params(6);   % angele beween v1_zrot and v2
            kappa1      = params(7);
            kappa2      = params(8);
                       
            % read intermediate variables
            cdtDiffPar = tmpAccessMemory.cdtDiffPar;
            cdtDiffPerp1 = tmpAccessMemory.cdtDiffPerp1;
            cdtDiffPerp2 = tmpAccessMemory.cdtDiffPerp2;
            cB = tmpAccessMemory.cB;
            gcB_gK1 = tmpAccessMemory.gcB_gK1;
            gcB_gK2 = tmpAccessMemory.gcB_gK2;
                
            % get coupled tensor jacobian
            cdt = tensor(); % Coupled Diffusion Tensor
            coupledTensorParams = [s0, cdtDiffPar, cdtDiffPerp1,...
                cdtDiffPerp2, theta, phi, alpha]';
            coupledTensorJac = cdt.jacobian(coupledTensorParams, schemefile);
            
            % initialise the jac with zeros
            jac = zeros(schemefile.measurementsNum(), obj.nParams);
                      
            % gradient wrt to s0
            jac(:,1) = coupledTensorJac(:,1);
            
            % gradient wrt dPar
            gf_gcdtDiffPar = coupledTensorJac(:,2);
            % gradient wrt dPerp1
            gf_gcdtDiffPerp1 = coupledTensorJac(:,3);
            % gradient wrt dPerp2
            gf_gcdtDiffPerp2 = coupledTensorJac(:,4);
            
            % gradient of dPar wrt diffPar
            gEcDiffPar_gDiffPar = -(gcB_gK1/cB);
            gEcDiffPerp1_gDiffPar = -(gcB_gK2/cB);
            gEcDiffPerp2_gDiffPar = 1 - gEcDiffPar_gDiffPar - gEcDiffPerp1_gDiffPar;
            
            jac(:,2) = gf_gcdtDiffPar.*gEcDiffPar_gDiffPar + gf_gcdtDiffPerp1.*gEcDiffPerp1_gDiffPar + gf_gcdtDiffPerp2.*gEcDiffPerp2_gDiffPar;
            
            % gradient of dPar wrt diffPerp
            gEcDiffPar_gDiffPerp = 1+(gcB_gK1/cB);
            gEcDiffPerp1_gDiffPerp = 1+(gcB_gK2/cB);
            gEcDiffPerp2_gDiffPerp = 2 - gEcDiffPar_gDiffPerp - gEcDiffPerp1_gDiffPerp;
            
            jac(:,3) = gf_gcdtDiffPar.*gEcDiffPar_gDiffPerp + gf_gcdtDiffPerp1.*gEcDiffPerp1_gDiffPerp + gf_gcdtDiffPerp2.*gEcDiffPerp2_gDiffPerp;
            
            % gradient wrt theta, phi, and alpha
            jac(:,4:6) = coupledTensorJac(:,5:7);
            
            % gradient wrt kappa1 and kappa2
            gcB_gK1gK2 = bingham.get_gcB_gk1gk2(0, kappa1, kappa2);
            gcB_gK1gK3 = bingham.get_gcB_gk1gk3(0, kappa1, kappa2);
            gcB_gK3 = bingham.get_gcB_gk3(0, kappa1, kappa2);
            gcB_gK2gK2 = bingham.get_gcB_gk2gk2(0, kappa1, kappa2);
            gcB_gK2gK3 = bingham.get_gcB_gk2gk3(0, kappa1, kappa2);
            %gcB_gK3gK3 = bingham.get_gcB_gk3gk3(0, kappa1, kappa2);
            
            gcdtDiffPar_gKappa1 = (diffPar-diffPerp)*((gcB_gK1/cB)*(gcB_gK2/cB)-gcB_gK1gK2/cB);
            gcdtDiffPar_gKappa2 = (diffPar-diffPerp)*((gcB_gK1/cB)*(gcB_gK3/cB)-gcB_gK1gK3/cB);
            
            gcdtDiffPerp1_gKappa1 =(diffPar-diffPerp)*((gcB_gK2/cB)^2-gcB_gK2gK2/cB);
            gcdtDiffPerp1_gKappa2 =(diffPar-diffPerp)*((gcB_gK2/cB)*(gcB_gK3/cB)-gcB_gK2gK3/cB); 
            
            %gcdtDiffPerp2_gKappa1 = (diffPar-diffPerp)*((gcB_gK2/cB)*(gcB_gK3/cB)-gcB_gK2gK3/cB);
            gcdtDiffPerp2_gKappa1 = - gcdtDiffPar_gKappa1 - gcdtDiffPerp1_gKappa1;
            %gcdtDiffPerp2_gKappa2 = (diffPar-diffPerp)*((gcB_gK3/cB)^2-gcB_gK3gK3/cB);
            gcdtDiffPerp2_gKappa2 = - gcdtDiffPar_gKappa2 - gcdtDiffPerp1_gKappa2;
            
            jac(:,7) = gf_gcdtDiffPar  .*gcdtDiffPar_gKappa1 + ...
                       gf_gcdtDiffPerp1.*gcdtDiffPerp1_gKappa1 + ...
                       gf_gcdtDiffPerp2.*gcdtDiffPerp2_gKappa1;
            
            
            jac(:,8) = gf_gcdtDiffPar  .*gcdtDiffPar_gKappa2 + ...
                       gf_gcdtDiffPerp1.*gcdtDiffPerp1_gKappa2 + ...
                       gf_gcdtDiffPerp2.*gcdtDiffPerp2_gKappa2;
         end % of jacobian  
        %\\
    end % of method (public) 
   %
end%of class zeppelinCBD
