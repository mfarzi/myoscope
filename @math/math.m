classdef math  
    % MATH 
    %
    %   math is a static class encapsulating constant variables and mehtods
    %   used for tensorial operations. This class is used with compartment
    %   models.
    %
    %
    %   methods (public):
    %       fit             - fit model parameters to measured signals
    %       fitMultiRun     - run the optimisation using different 
    %                         initial points selected randomly.
    %       synthesize      - synthesize attenuation signals
    %       getCost         - return the Root Mean Square (RMS) error of
    %                         the fitted model
    %       getParams       - return the model parameters in a column
    %                         vector
    %       getParamsNum    - return the number of model parameters
    %       getFixedParams  - return a logical column vector stating if a
    %                         parameter is fixed (true) or not (false).
    %
    %
    %   See also: compartment
    %
    % Mohsen Farzi
    % Email: m.farzi@leeds.ac.uk
    
    properties (Constant=true)
        % the gyromagnetic Ratio (copied from camino source codes)
        GAMMA = 2.6751525e8; 
        
        % 60 first roots from the equation J'1(x)=0
        % J'1 is the derivative of the first order bessel function
        % (copied from camino source codes)
        Jp1ROOTS = ...
         [1.84118307861360, 5.33144196877749, 8.53631578218074, ...
          11.7060038949077, 14.8635881488839, 18.0155278304879, ...
          21.1643671187891, 24.3113254834588, 27.4570501848623, ...
          30.6019229722078, 33.7461812269726, 36.8899866873805, ...
          40.0334439409610, 43.1766274212415, 46.3195966792621, ...
          49.4623908440429, 52.6050411092602, 55.7475709551533, ...
          58.8900018651876, 62.0323477967829, 65.1746202084584, ...
          68.3168306640438, 71.4589869258787, 74.6010956133729, ...
          77.7431620631416, 80.8851921057280, 84.0271895462953, ...
          87.1691575709855, 90.3110993488875, 93.4530179063458, ...
          96.5949155953313, 99.7367932203820, 102.878653768715, ...
          106.020498619541, 109.162329055405, 112.304145672561, ...
          115.445950418834, 118.587744574512, 121.729527118091, ...
          124.871300497614, 128.013065217171, 131.154821965250, ...
          134.296570328107, 137.438311926144, 140.580047659913, ...
          143.721775748727, 146.863498476739, 150.005215971725, ...
          153.146928691331, 156.288635801966, 159.430338769213, ...
          162.572038308643, 165.713732347338, 168.855423073845, ...
          171.997111729391, 175.138794734935, 178.280475036977, ...
          181.422152668422, 184.563828222242, 187.705499575101]';
    end
    
    methods (Static)
        leb_tmp = getLebedevSphere(degree);
        U = getUnitFrame(theta, phi, alpha);
        u1 = get_u1(theta, phi, ~);
        gU1_gTheta = get_gU1_gTheta(theta, phi, ~)
        gU1_gPhi = get_gU1_gPhi(theta, phi, ~)
        [theta, phi, alpha] = getOrientationAngles(U);
        [gU1, gU2, gU3] = getOrientationJacobian(theta, phi, alpha);
        state = isOrthonormal(U);
        function rotateAxis(obj)
            % rotateAxis remove the ambiguity in estimation of parameters
            % theta, phi, and alpha by rotating the cooriante system of the
            % eigen vectors appropriately.
            % 
            % The orthonormal coordiante system V = [v1, v2, v3] should 
            % rotate such that the eigne vectors v1, v2, and v3 are 
            % associated with the eigen values in a descending order.
            %
            % More over, the direction of v1 or -v1 is selected such that
            % the angle theta is always in range [0, pi/2].
            % 
            % Finally, the direction v2 or -v2 is selected such that the
            % angle alpha is always between [0, phi].
            %
            % see also: getUnitFrame, getEulerAngles
            
            % read params
            p = obj.params;
            diffPar     = p(2);   % diffusivity [s/m^2]
            diffPerp1   = p(3);   % diffusivity [s/m^2]
            diffPerp2   = p(4);   % diffusivity [s/m^2]
            theta       = p(5);   % elevation angle
            phi         = p(6);   % azimuth angel
            alpha       = p(7);   % angele beween v1_zrot and v2
            
            V = math.getUnitFrame(theta, phi, alpha);
            
            % sort the eigen values in descending order
            eigVal = [diffPar; diffPerp1; diffPerp2];
            [eigVal, sortID] = sort(eigVal, 'descend');
            
            diffPar   = eigVal(1);
            diffPerp1 = eigVal(2);
            diffPerp2 = eigVal(3);
            
            V = V(:,sortID);
            % make sure V is still orthogonal following the permutation
            if uint8(norm(cross(V(:,1), V(:,2))-V(:,3))) ~= 0
                V(:,3) = -V(:,3);
            end
            [theta, phi, alpha] = math.getEulerAngles(V);
            
            % make sure theta < pi/2
            if obj.theta>pi/2
                theta = pi - theta;
                phi = pi + phi;
                alpha = -alpha;
            end

            phi = mod(phi, 2*pi);
            
            % make sure alpha is in rage [0, pi]
            alpha = mod(alpha, pi);
            
            obj.params(2:end) = [diffPar; diffPerp1; diffPerp2; theta; ...
                                 phi; alpha];
        end % of rotateAxis
        
        function adc = getADC(obj)
            % compute Apparate Diffusion Coefficient
            adc = (obj.diffPar + obj.diffPerp1 + obj.diffPerp2)/3;
        end
        
        function fa = getFA(obj)
            % Fractional anisotropy
            lambda = [obj.diffPar, obj.diffPerp1, obj.diffPerp2];
            lambda_hat = mean(lambda);
            fa = sqrt(3/2)*norm(lambda-lambda_hat)/norm(lambda);
        end
                    
        function DT = getDT(obj)
            V = obj.getEigenVec();
            D = diag([obj.diffPar; obj.diffPerp1; obj.diffPerp2]);
            DT = V*D*V';
        end
        
        function V = getEigenVec(obj)
            V = getUnitFrame(obj.theta, obj.phi, obj.alpha);
        end
    end%of method
end