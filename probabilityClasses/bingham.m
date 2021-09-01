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
        
        
        
        
        function [gcB_gK1, gcB_gK2] = cBinghamGrad(k1, k2)
            gcB_gK1 = integral2(@(theta,phi) ...
                                exp(k1*sin(theta).^2.*cos(phi).^2  + ...
                                k2*sin(theta).^2.*sin(phi).^2).* ...
                                sin(theta).^3.*cos(phi).^2, ...
                                0,pi,0,2*pi)/(4*pi);
                            
            gcB_gK2 = integral2(@(theta,phi) ...
                                exp(k1*sin(theta).^2.*cos(phi).^2  + ...
                                k2*sin(theta).^2.*sin(phi).^2).* ...
                                sin(theta).^3.*sin(phi).^2, ...
                                0,pi,0,2*pi)/(4*pi);                
        end
        
        function g_kappa = hypergeomGrad(a, b, kappa)
            % M(a, b, kappa)= (gamma(b)/(gamma(b-a)*gamma(a))) ... 
            %           * integral(@(t) exp(kappa.*t).*t.^(a-1).*(1-t).^(b-a-1),0,1);
            g_kappa = (gamma(b)/(gamma(b-a)*gamma(a))) ...
                    * integral(@(t) exp(kappa.*t).*t.^(a).*(1-t).^(b-a-1),0,1);
        end
        
        function p = binghamPdf(n, u1, u2, kappa1, kappa2, cB)
            p = exp(kappa1*(u1'*n).^2 + kappa2*(u2'*n).^2)/cB;
        end
        
    end
    
    methods (Static)
        function probs = pdf(params, orientations)
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
        
        function [samples, probs] = samplepdf(params, n)
            % samplepdf(params, n) return n points on the unit sphere with
            % their corresponding probability weight normalised to one.
            %
            %   input arguments:
            %            params: A numeric column vector of size 6
            %                   [theta; phi; alpha; kappa1; kappa2; kappa3]
            %                 n: Positive integer specifying number of 
            %                    points on the unit sphere.
            %                    Allowed values: {   6,  14,  26,  38,  50,
            %                                       74,  86, 110, 146, 170,
            %                                      194, 230, 266, 302, 350,
            %                                      434, 590, 770, 974,1202,
            %                                     1454,1730,2030,2354,2702,
            %                                     3074,3470,3890,4334,4802,
            %                                     5294,5810}
            %                    see math.getLebedevSphere
            %
            %  output arguments:
            %           samples: A struct with fields x, y, z, and w
            %                    returning n selected orientations on the 
            %                    unit sphere and the normalising weights to
            %                    the unit sphere area (4*pi)
            %             probs: A numeric column vector of size n
            %                    returning probabilities for each selected 
            %                    orientation drawn from the Bingham
            %                    distribution.  
            %
            % get samples
            samples = math.getLebedevSphere(n);
            orientations = [samples.x, samples.y, samples.z];
            
            % get probs
            probs = bingham.pdf(params, orientations);
            probs = probs(:).*samples.w/(4*pi);
        end
    end%of methods (Static)
end