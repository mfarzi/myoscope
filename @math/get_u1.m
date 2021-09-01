function u1 = get_u1(theta, phi, ~)
    % get_u1 return the primary eigenvector from an orthonormal coordinate
    % system U = [u1, u2, u3] represented by angles theta, phi, and alpha.
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
    
    assert(isrow(theta), 'MATLaB:math:get_u1',...
        'theta must be a row vector');
    
    assert(isrow(phi), 'MATLaB:math:get_u1',...
        'phi must be a row vector');
    
    assert(length(theta)==length(phi), 'MATLaB:math:get_u1',...
        'theta and phi must be of the same size.');
    
    u1 = [sin(theta).*cos(phi); ...
          sin(theta).*sin(phi); ...
          cos(theta)];