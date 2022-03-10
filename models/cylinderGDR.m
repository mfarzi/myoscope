classdef cylinderGDR < compartment
    % CYLINDERGDR
    % (Cylinder with Gamma Distribution Radii)
    %
    %   An CYLINDERGDR object is a basic COMPARTMENT object represenitng 
    %   particles diffusion in a cylinder with various radii drawn from a 
    %   Gamma distribution. This class implements inherited abstract 
    %   methods and properties. This class should be used as input to
    %   MULTICOMPARTMENT class.
    %
    %   Model Parameters:
    %       name              - class name
    %       s0                - Normalised b0-signal 
    %       diff              - diffusivity along the cylinder axis [m^2/s]
    %       kappa             - shape parameter for gamma distribution
    %       nu                - scale parameter for gamma distribution
    %       theta             - elevation angle; the angle between the
    %                           cylinder axis (n1) and the z-axis. [radian]
    %       phi               - azimuth angle; the angle between the x-axis
    %                           and the n1 projection onto the xy-plane. 
    %                           [radian]
    %
    %   For mathematical background see
    %       Y. Assaf, R. Z. Freidlin, G. K. Rohde, and P. J. Basser,
    %       New modeling and experimental framework to characterize 
    %       hindered and restricted water diffusion in brain white 
    %       matter, Magn. Reson. Med., vol. 52, no. 5, pp. 965?978, 2004.
    %
    %   See also: multicompartment, cylinder, cylinderECS, cylinderBDA,
    %             stick
    %
    % Mohsen Farzi
    % Email: m.farzi@leeds.ac.uk
    
    properties (SetAccess='protected')
        name = 'cylinderGDR';    
        nCompartments = 1;     % number of basic COMPARTMETNT objects
        nParams = 6;
        nHyperparams = 1;
    end
    
    properties (Access = 'protected')
        comp = [];             % List of COMPARTMENT object 
        links;
        hyperparams = [];
        hyperparamsName = {};
    end
    
    methods
        function obj = cylinderGDR()
            %CYLINDERGDR Construct Function.
            %
            %   cylinderGDR() constructs a single-compartment model
            
            % initialise the linking functions to map constrained model
            % parameters to unconstrained optimisation variables.
            s0      = linker('s0'       , 'bounded', 0    , 1);
            diffPar = linker('diffPar'  , 'bounded', 1e-10, 3e-9);
            kappa   = linker('kappa'    , 'bounded', 1, 100);
            nu      = linker('nu'       , 'bounded', 0.1e-6, 10e-6);
            theta   = linker('theta');
            phi     = linker('phi');
            obj.links   = [s0; diffPar; kappa; nu; theta; phi];
            
            % set hparams and hparams Names
            obj.hyperparams = 30;
            obj.hyperparamsName = {'nBinsR'};
        end
        
        function s = synthesize(obj, params, schemefile)
            % synthesize(params, scheme, hyperparams) return DW-MR signals
            % from a cylinder with Gamma Distributed Radii. The signal 
            % is modeled as the average of the signals from a cylinder with 
            % radius r where r is drawn from Gamma(kappa, nu).
            %       Input arguments:
            %                params: Numerical column vector of all model
            %                        parameters
            %                scheme: Diffusion scheme
            %           hyperparams: nBinsR=50; number of points to compute
            %                        the average signal
            %
            %      output arguments:
            %                     s: Numerical column vector of synthesied
            %                        signal
            %  
            
            % validate inputs
            validateattributes(params, {'double'},...
                {'column', 'nrows', obj.nParams},...
                'cylinderGDR.synthesize', 'params');
            
            assert(isa(schemefile, 'scheme'), ...
                'MATLAB:tensor:invalidInputArgument',...
                'Scheme file should be of type scheme.');
            
            assert(strcmp(schemefile.type, 'stejskal-tanner'), ...
                'MATLAB:tensor:invalidInputArgument',...
                'Scheme file should be of type stejskal-tanner.');
            
             % read params into individual model parameters
            s0      = params(1);   % b0 signal
            diffPar = params(2);   % diffusivity [s/m^2]
            kappa   = params(3);   % shape parameter for Gamma
            nu      = params(4);   % scale parameter for Gamma
            theta   = params(5);   % elevation angel [radian]
            phi     = params(6);   % azimuth angle [radian]
            
            % read hyper-parameters
            nBinsR = obj.hyperparams(1);
            
            % discritize the radii
            [probs, radii] = cylinderGDR.discritiseGammaDistribution(kappa, nu, nBinsR);
            
            P = ones(5, nBinsR);
            P(1,:) = s0;
            P(2,:) = diffPar;
            P(3,:) = radii';
            P(4,:) = theta;
            P(5,:) = phi;
            
            tmpCylinder = cylinder();
            tmp = arrayfun(@(n) tmpCylinder.synthesize(P(:,n), schemefile), 1:nBinsR, 'UniformOutput', false);
            sig = cell2mat(tmp);
            
            % compute the average wrt R
            s = sig*probs;
        end
        
        function jac = jacobian(obj, params, schemefile)
            % jacobian(params, scheme, hyperparams) return the gradient of 
            % signal wrt to model parameters.
            %
            %       Input arguments:
            %                params: Numerical column vector of all model
            %                        parameters
            %                scheme: Diffusion scheme
            %           hyperparams: Not used in this class. Reserved 
            %                        variable for non-optimisable model
            %                        parameters.
            %
            %      output arguments:
            %                   jac: Numerical matrix of size M x nParams.
             
            % validate inputs
            validateattributes(params, {'double'},...
                {'column', 'nrows', obj.nParams},...
                'cylinderGDR.jacobian', 'params');
            
            assert(isa(schemefile, 'scheme'), ...
                'MATLAB:cylinderGDR:invalidInputArgument',...
                'Scheme file should be of type scheme.');
            
            assert(strcmp(schemefile.type, 'stejskal-tanner'), ...
                'MATLAB:cylinderGDR:invalidInputArgument',...
                'Scheme file should be of type stejskal-tanner.');
            
            % read params into individual model parameters
            s0      = params(1);   % b0 signal
            diffPar = params(2);   % diffusivity [s/m^2]
            kappa   = params(3);   % shape parameter for Gamma
            nu      = params(4);   % scale parameter for Gamma
            theta   = params(5);   % elevation angel [radian]
            phi     = params(6);   % azimuth angle [radian]
            
            % read hyper-parameters
            nBinsR = obj.hyperparams(1);
            
            % discritize the radii
            [probs, radii] = cylinderGDR.discritiseGammaDistribution(kappa, nu, nBinsR);
            
            
            P = ones(5, nBinsR);
            P(1,:) = s0;
            P(2,:) = diffPar;
            P(3,:) = radii;
            P(4,:) = theta;
            P(5,:) = phi;
            
            tmpCylinder = cylinder();
            [tmp1, tmp2] = arrayfun(@(n) tmpCylinder.jacobian(P(:,n), schemefile), 1:nBinsR, 'UniformOutput', false);
            jac_r = cell2mat(tmp1);
            sig = cell2mat(tmp2);
            N = nBinsR*tmpCylinder.nParams;
            
            gSgS0      = jac_r(:,1:tmpCylinder.nParams:N)*probs;
            gSgDiffPar = jac_r(:,2:tmpCylinder.nParams:N)*probs;
            gSgKappa   = sig*((log(radii)-log(nu)-psi(kappa)).*probs);
            gSgNu      = sig*(probs.*(radii/nu^2 - kappa/nu));
            gSgTheta   = jac_r(:,4:tmpCylinder.nParams:N)*probs;
            gSgPhi     = jac_r(:,5:tmpCylinder.nParams:N)*probs;
            
            
            jac = [gSgS0, gSgDiffPar, gSgKappa, gSgNu, gSgTheta, gSgPhi];
        end
        %\\
    end % of methods (public)
    
    methods (Static, Access = 'protected')
         function [probs, radii] = discritiseGammaDistribution(kappa, nu, nBinsR)
            % discritiseGammaDistribution returns a pair of probabilities 
            % and corresponding radii for numerical computation of 
            % intergrals.
            
            if nBinsR == 1
                radii = nu * kappa;
                probs = 1;
                return;
            end
            
            rMin = icdf('gamma', 1/nBinsR, kappa, nu);
            %rMin = 0;
            
            rMax = icdf('gamma', 1-1/nBinsR, kappa, nu);
            %rMax = 30e-6;
            
            r = linspace(rMin, rMax, nBinsR+1)';
            
            % compute probabilitis in each bin
            probs = diff(cdf('gamma', r, kappa, nu)); 
            
            % normalise probabilities
            probs = probs/sum(probs);
            
            % compute average radius in each bin
            radii = (r(1:nBinsR) + r(2:nBinsR+1))*0.5;
            %probs = pdf('gamma', radii, kappa, nu);
         end  
    end
end % of cylinder class

