function [theta, phi, alpha] = getEulerAngles(U)
    % getEulerAngles return the required rotation angles to transform the
    % standard orthonormal coordinate system to yield the given coordinate 
    % system U = [u1, u2, u3].
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
    % see also: getUnitFrame
    %
    % Mohsen Farzi
    % m.farzi@leeds.ac.uk
    
    if ~(all(size(U) == [3, 3]) && isnumeric(U))
        error('The input matrix should be of size 3x3.\n');
    end
    
    u1 = U(:,1); u2 = U(:,2); u3 = U(:, 3);
    
    % check if the input frame is normal
    if not(uint8(norm(u1))==1 && uint8(norm(u2))==1 && uint8(norm(u3))==1)
        error('Columns of the input matrix are not normal.\n');
    end
    
    % check if the input frame is orthogonal
    if uint8(norm(cross(u1,u2)-u3))==2
        u3 = -u3;
    elseif uint8(norm(cross(u1,u2)-u3)) ~= 0
        error(['Columns of the input matrix are not orthogonal.\n',...
               'norm(u1Xu2-u3)=%f\n'], norm(cross(u1, u2)-u3));
    end
    
    % esitmate theta, phi, alpha
    theta = acos(u1(3));
    sinTheta = sqrt(1-u1(3)^2);

    if sinTheta~= 0
        phi = atan2(u1(2)/sinTheta, u1(1)/sinTheta);
        alpha = atan2(u3(3)/sinTheta, -u2(3)/sinTheta);
    else
        % if theta = 0, then one degree of freedom is enough. so, phi = 0.
        phi = 0;
        alpha = atan2(u2(2), u2(1));
    end
    
end


