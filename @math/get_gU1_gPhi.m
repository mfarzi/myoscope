function gU1_gPhi = get_gU1_gPhi(theta, phi, ~)
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
    
    assert(isrow(theta), 'MATLaB:math:get_gU1_gTheta',...
        'theta must be a row vector');
    
    assert(isrow(phi), 'MATLaB:math:get_gU1_gTheta',...
        'phi must be a row vector');
    
    assert(length(theta)==length(phi), 'MATLaB:math:get_gU1_gTheta',...
        'theta and phi must be of the same size.');
    
    %% jacobian of u1 wrt theta, phi, and alpha
    %     u1 = [sin(theta)*cos(phi); ...
    %           sin(theta)*sin(phi); ...
    %           cos(theta)];
    z = zeros(1,length(theta));
    gU1_gPhi   = [-sin(phi).*sin(theta); ...
                   cos(phi).*sin(theta); ...
                   z]; 