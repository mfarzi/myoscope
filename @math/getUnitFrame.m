function U = getUnitFrame(theta, phi, alpha)
    % getUnitFrame return an orthonormal coordinate system U = [u1, u2, u3]
    % from Euler angles theta, phi, and alpha.
    %
    % To define the input angles, note that the vector n1 is defined as
    % the vector u1 rotated by pi/2 around the z-axis. The vector p1 is the
    % projection of u1 on the xy-plane.
    %
    %       theta:         the angle between z-axis and u1. This angle is
    %                      in the range [0, pi].
    %
    %         phi:         the angle between p1 and x-axis. This angle is
    %                      in the range [-pi, pi] using right-hand rule.
    %
    %       alpha:         the angle between n1 and u2. This angle is
    %                      between [-pi, pi] using right-hand rule. 
    %
    % NOTE: THE EULER ANGLES DEFINED HERE ARE DIFFERENT FROM THE STANDARD
    % DEFINITION. HOWEVER, THEY ARE CALLED EULER ANGLES AS THE UNDERLYING
    % CONCEPT IS SIMILAR.
    %
    % see also: getEulerAngles
    %
    % Mohsen Farzi
    % m.farzi@leeds.ac.uk

    u1 = [sin(theta)*cos(phi); ...
          sin(theta)*sin(phi); ...
          cos(theta)];
                
    % rotate u1 by pi/2 around the xy-plane
    % n1 = [ sin(theta+pi/2)*cos(phi); ...
    %        sin(theta+pi/2)*sin(phi); ...
    %        cos(theta+pi/2)];
    % rotate around v1 by alpha
    % u2 = rodrigues_rot(n1, u1, alpha);
    u2 = [ cos(theta)*cos(alpha)*cos(phi) - sin(alpha)*sin(phi); ...
           cos(theta)*cos(alpha)*sin(phi) + sin(alpha)*cos(phi); ...
          -sin(theta)*cos(alpha)];
              
    % use cross vector product to compute v3
    % u3 = cross(u1, u2);
    u3 = [-cos(theta)*sin(alpha)*cos(phi) - cos(alpha)*sin(phi); ...
          -cos(theta)*sin(alpha)*sin(phi) + cos(alpha)*cos(phi); ...
           sin(theta)*sin(alpha)];
       
    U = [u1, u2, u3];
end