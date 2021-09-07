classdef bingham
    % BINGHAM
    %   A bingham object compute the bingham distribution.
    %
    %   For mathematical background see
    %   C. Bingham, An antipodally symmetric distribution on the sphere,
    %   The Annals of Statistics, vol. 2, no. 6, pp. 1201-1225, 1974
    %
    % Mohsen Farzi
    % Email: m.farzi@leeds.ac.uk
    
    properties (Access = 'protected')
        links;
    end
    
    
    methods
        function g_kappa = hypergeomGrad(a, b, kappa)
            % M(a, b, kappa)= (gamma(b)/(gamma(b-a)*gamma(a))) ... 
            %           * integral(@(t) exp(kappa.*t).*t.^(a-1).*(1-t).^(b-a-1),0,1);
            g_kappa = (gamma(b)/(gamma(b-a)*gamma(a))) ...
                    * integral(@(t) exp(kappa.*t).*t.^(a).*(1-t).^(b-a-1),0,1);
        end
    end
    
    methods (Static)
        function probs = pdf(params, orientations, weights)
            % pdf(params, orientations) return probability density function
            % at given orientations
            %
            %   input arguments:
            %            params: A numeric column vector of size 6
            %                   [theta; phi; alpha; kappa1; kappa2; kappa3]
            %      orientations: A numeric matrix of size [npts x 3]. Each
            %                    column represents the x, y, and z 
            %                    coordinates respectively.
            %
            %  output arguments:
            %             probs: A numeric column vector of size npts
            %                    returning probabilities for each selected 
            %                    orientation drawn from the Bingham
            %                    distribution. 
            %
            
            % read Bingham parameters
            theta  = params(1);
            phi    = params(2);
            alpha  = params(3);
            kappa1 = params(4);
            kappa2 = params(5);
            kappa3 = params(6);
            
            U = math.getUnitFrame(theta, phi, alpha);
            x = orientations*U;
            probs = exp(-kappa1*x(:,1).^2 - kappa2*x(:,2).^2 - kappa3*x(:,3).^2);
            
            c = bingham.getC(kappa1, kappa2, kappa3);
            probs = probs/c;
            probs = probs.*weights;
        end
        
        function c = getC(kappa1, kappa2, kappa3)
            % returns the hypergeometric function of the first kind
            % c = 1F1(1/2,3/2,B) = integral_S^2 exp(n'Bn) dn
            % It can be shown that B = U*Z*U', then c is independent of U.
                 
            c = integral2(@(theta,phi) ...
                exp(-kappa1*sin(theta).^2.*cos(phi).^2  - ...
                     kappa2*sin(theta).^2.*sin(phi).^2  - ...
                     kappa3*cos(theta).^2).* sin(theta),0,pi,0,2*pi)/(4*pi);      
            
        end
        
        function jac = jacobian(params, orientations, weights)
            % joacobian(params, orientations) return the gradient of P_B
            % wrt Bingham parameters
            %
            %   input arguments:
            %            params: A numeric column vector of size 6
            %                   [theta; phi; alpha; kappa1; kappa2; kappa3]
            %      orientations: A numeric matrix of size [npts x 3]. Each
            %                    column represents the x, y, and z 
            %                    coordinates respectively.
            %
            %  output arguments:
            %               jac: A numerical matrix of size npts x 6
            %
            
            % compute probabilities at each orientation
            probs = bingham.pdf(params, orientations, weights);
            
            % read Bingham parameters
            theta  = params(1);
            phi    = params(2);
            alpha  = params(3);
            kappa1 = params(4);
            kappa2 = params(5);
            kappa3 = params(6);
            
            npts = size(orientations, 1);
            jac = zeros(npts, 6);
            
            U = math.getUnitFrame(theta, phi, alpha);
            UdotO = orientations*U;
            
            gPb_gu1 = -2*kappa1*UdotO(:,1).*probs.*orientations;
            gPb_gu2 = -2*kappa2*UdotO(:,2).*probs.*orientations;
            gPb_gu3 = -2*kappa3*UdotO(:,3).*probs.*orientations;
            
            [gU1, gU2, gU3] = math.getOrientationJacobian(theta, phi, alpha);
            
            % gPb_gTheta
            jac(:,1) = gPb_gu1*gU1(:,1) + gPb_gu2*gU2(:,1) + gPb_gu3*gU3(:,1);   
            % gPb_gPhi
            jac(:,2) = gPb_gu1*gU1(:,2) + gPb_gu2*gU2(:,2) + gPb_gu3*gU3(:,2); 
            % gPb_gAlpha
            jac(:,3) = gPb_gu1*gU1(:,3) + gPb_gu2*gU2(:,3) + gPb_gu3*gU3(:,3); 
            
            
            cB = bingham.getC(kappa1, kappa2, kappa3);
            
            gcB_gK1 = bingham.get_gcB_gk1(kappa1, kappa2, kappa3); 
                       
            gcB_gK2 = bingham.get_gcB_gk2(kappa1, kappa2, kappa3);
                     
            gcB_gK3 = bingham.get_gcB_gk3(kappa1, kappa2, kappa3);
            %gPb_gK1
            jac(:,4) = -probs.*(UdotO(:,1).^2+gcB_gK1/cB);
            jac(:,5) = -probs.*(UdotO(:,2).^2+gcB_gK2/cB);
            jac(:,6) = -probs.*(UdotO(:,3).^2+gcB_gK3/cB);
        end
        
        function g = get_gcB_gk1(kappa1, kappa2, kappa3)
            g = integral2(@(theta,phi) ...
                -exp(-kappa1*sin(theta).^2.*cos(phi).^2  - ...
                      kappa2*sin(theta).^2.*sin(phi).^2  - ...
                      kappa3*cos(theta).^2).* sin(theta).^3.*cos(phi).^2,0,pi,0,2*pi)/(4*pi);
        end
        
        function g = get_gcB_gk1gk2(kappa1, kappa2, kappa3)
            g = integral2(@(theta,phi) ...
                exp(-kappa1*sin(theta).^2.*cos(phi).^2  - ...
                      kappa2*sin(theta).^2.*sin(phi).^2  - ...
                      kappa3*cos(theta).^2).* sin(theta).^5.*cos(phi).^2.*sin(phi).^2,0,pi,0,2*pi)/(4*pi);
        end
        
        function g = get_gcB_gk1gk3(kappa1, kappa2, kappa3)
            g = integral2(@(theta,phi) ...
                exp(-kappa1*sin(theta).^2.*cos(phi).^2  - ...
                      kappa2*sin(theta).^2.*sin(phi).^2  - ...
                      kappa3*cos(theta).^2).* sin(theta).^3.*cos(phi).^2.*cos(theta).^2,0,pi,0,2*pi)/(4*pi);
        end
        
        function g = get_gcB_gk2(kappa1, kappa2, kappa3)
            g = integral2(@(theta,phi) ...
                -exp(-kappa1*sin(theta).^2.*cos(phi).^2  - ...
                      kappa2*sin(theta).^2.*sin(phi).^2  - ...
                      kappa3*cos(theta).^2).* sin(theta).^3.*sin(phi).^2,0,pi,0,2*pi)/(4*pi);
        end
        
        function g = get_gcB_gk2gk2(kappa1, kappa2, kappa3)
            g = integral2(@(theta,phi) ...
                exp(-kappa1*sin(theta).^2.*cos(phi).^2  - ...
                      kappa2*sin(theta).^2.*sin(phi).^2  - ...
                      kappa3*cos(theta).^2).* sin(theta).^5.*sin(phi).^4,0,pi,0,2*pi)/(4*pi);
        end
        
        function g = get_gcB_gk2gk3(kappa1, kappa2, kappa3)
            g = integral2(@(theta,phi) ...
                exp(-kappa1*sin(theta).^2.*cos(phi).^2  - ...
                      kappa2*sin(theta).^2.*sin(phi).^2  - ...
                      kappa3*cos(theta).^2).* sin(theta).^3.*sin(phi).^2.*cos(theta).^2,0,pi,0,2*pi)/(4*pi);
        end
        
        function g = get_gcB_gk3(kappa1, kappa2, kappa3)
            g = integral2(@(theta,phi) ...
                -exp(-kappa1*sin(theta).^2.*cos(phi).^2  - ...
                      kappa2*sin(theta).^2.*sin(phi).^2  - ...
                      kappa3*cos(theta).^2).* sin(theta).*cos(theta).^2,0,pi,0,2*pi)/(4*pi);
        end
        
        function g = get_gcB_gk3gk3(kappa1, kappa2, kappa3)
            g = integral2(@(theta,phi) ...
                exp(-kappa1*sin(theta).^2.*cos(phi).^2  - ...
                      kappa2*sin(theta).^2.*sin(phi).^2  - ...
                      kappa3*cos(theta).^2).* sin(theta).*cos(theta).^4,0,pi,0,2*pi)/(4*pi);
        end
        
        function B = getB(params)
            % return the input params in the format of matrix B
            theta  = params(1);
            phi    = params(2);
            alpha  = params(3);
            kappa1 = params(4);
            kappa2 = params(5);
            kappa3 = params(6);
            
            U = math.getUnitFrame(theta, phi, alpha);
            B = U*diag([kappa1, kappa2, kappa3])*U';
        end
    end%of methods (Static)
end