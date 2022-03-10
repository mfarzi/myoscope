function [theta, phi, alpha] = simplifyOrientationAngles(theta, phi, alpha)
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
    N = length(theta);
    assert(length(phi)==N && length(alpha)==N,...
        'MATLAB:math:invalidInputArgs', ...
        'Input vectors must be of the same lenght.');
    
    theta = mod(theta, 2*pi);
    phi = mod(phi, 2*pi);
    alpha = mod(alpha, 2*pi);
    
    % theta should be periodic with pi
    idx = theta>pi;
    theta(idx) = theta(idx)-pi;
    alpha(idx) = 2*pi-alpha(idx);
    alpha = mod(alpha, pi);
    
    % make sure theta < pi/2 (for uniqueness of solution)
    idx = theta>pi/2;
    theta(idx) = pi-theta(idx);
    phi(idx) = pi + phi(idx);
    phi = mod(phi, 2*pi);
    phi(phi>pi) = phi(phi>pi)-2*pi;
    alpha(idx) = pi-alpha(idx);
    
    
end