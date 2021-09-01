function [theta, phi, alpha] = getOrientationAngles(U)
    % getOrientationAngles represent an input orthonormal coordinate system
    % by means of three rotation angels theta, phi, and alpha. 
    %
    % To define the input angles, note that the vector n1 is defined as
    % the vector u1 rotated by pi/2 around the z-axis. The vector p1 is the
    % projection of u1 on the xy-plane.
    %       Input Arguments:
    %                     U: An orthonormal basis matrix of 3x3.
    %                        U = [u1, u2, u3]construct the frame.
    %      Output Arguments:
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
    %
    % see also: getUnitFrame, getTensorialAngles
    %
    % Mohsen Farzi
    % m.farzi@leeds.ac.uk
    
    % check input size and type
    validateattributes(U, {'numeric'}, {'size', [3,3]},...
        'getOrientationAngles', 'U');
    
    % check if input matrix is orthonormal
    assert(math.isOrthonormal(U), 'MATLAB:math:getOrientationAngles',...
        "Input matrix U must be orthonormal,i.e U'U=eye(3)");
    
    % esitmate theta, phi, alpha
    u1 = U(:,1); u2 = U(:,2); u3 = U(:, 3);
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




