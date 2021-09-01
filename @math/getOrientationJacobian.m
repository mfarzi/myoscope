function [gU1, gU2, gU3] = getOrientationJacobian(theta, phi, alpha)
    % getOrientationJacobian return the jacobian matrix [3x3] for each 
    % unit vector with respect to orientation angles theta, phi, and alpha.
    %
    %       input arguments:
    %                 theta: The angle between z-axis and u1. 
    %                        Range [0, pi].
    %
    %                   phi: The angle between x-axis and p1 (u1 projection 
    %                        on xy-plane). 
    %                        Range: [-pi, pi] radian
    %
    %                 alpha: The angle between n1 (u1 rotated on zp1-plane 
    %                        by pi/2) and u2. 
    %                        Range: [-pi, pi] radian
    %      output arguments:
    %                   gu1: 3x3 Jacobian matrix for u1. Each row is the
    %                        gradient of u1(i) with respect to theta, phi,
    %                        and alpha.
    %                   gu2: 3x3 Jacobian matrix for u2. Each row is the
    %                        gradient of u2(i) with respect to theta, phi,
    %                        and alpha.
    %                   gu3: 3x3 Jacobian matrix for u3. Each row is the
    %                        gradient of u3(i) with respect to theta, phi,
    %                        and alpha.
    %                        
    % see also: getOrientationAngles, getUnitFrame, getTensorialAngles
    %
    % Mohsen Farzi
    % m.farzi@leeds.ac.uk
    
    %% jacobian of u1 wrt theta, phi, and alpha
    %     u1 = [sin(theta)*cos(phi); ...
    %           sin(theta)*sin(phi); ...
    %           cos(theta)];
    gU1_gTheta = [ cos(phi)*cos(theta); ...
                   sin(phi)*cos(theta); ...
                  -sin(theta)];
    gU1_gPhi   = [-sin(phi)*sin(theta); ...
                   cos(phi)*sin(theta); ...
                   0];           
    gU1_gAlpha = zeros(3, 1);
    gU1 = [gU1_gTheta, gU1_gPhi, gU1_gAlpha];
 
    %% jacobian of u2 wrt theta, phi, and alpha
    % u2 = [ cos(theta)*cos(alpha)*cos(phi) - sin(alpha)*sin(phi); ...
    %        cos(theta)*cos(alpha)*sin(phi) + sin(alpha)*cos(phi); ...
    %       -sin(theta)*cos(alpha)];
    gU2_gTheta = [-cos(alpha)*cos(phi)*sin(theta); ...
                  -cos(alpha)*sin(phi)*sin(theta); ...
                  -cos(alpha)*cos(theta)];
    gU2_gPhi   = [-sin(phi)*cos(alpha)*cos(theta)-cos(phi)*sin(alpha); ...          
                   cos(phi)*cos(alpha)*cos(theta)-sin(phi)*sin(alpha); ...
                   0];          
    gU2_gAlpha = [-sin(alpha)*cos(theta)*cos(phi)-cos(alpha)*sin(phi);...
                  -sin(alpha)*cos(theta)*sin(phi)+cos(alpha)*cos(phi);...
                   sin(alpha)*sin(theta)];           
    gU2 = [gU2_gTheta, gU2_gPhi, gU2_gAlpha];
    
    %% jacobian of u2 wrt theta, phi, and alpha
    % u3 = [-cos(theta)*sin(alpha)*cos(phi) - cos(alpha)*sin(phi); ...
    %       -cos(theta)*sin(alpha)*sin(phi) + cos(alpha)*cos(phi); ...
    %        sin(theta)*sin(alpha)];
    gU3_gTheta = [ sin(theta)*cos(phi)*sin(alpha);...
                   sin(theta)*sin(phi)*sin(alpha);...
                   cos(theta)*sin(alpha)];
    gU3_gPhi   = [ sin(phi)*cos(theta)*sin(alpha)-cos(phi)*cos(alpha);...
                  -cos(phi)*cos(theta)*sin(alpha)-sin(phi)*cos(alpha);...
                   0];    
    gU3_gAlpha = [-cos(alpha)*cos(theta)*cos(phi)+sin(alpha)*sin(phi);...
                  -cos(alpha)*cos(theta)*sin(phi)-sin(alpha)*cos(phi);...
                   cos(alpha)*sin(theta)];   
    gU3 = [gU3_gTheta, gU3_gPhi, gU3_gAlpha];
end