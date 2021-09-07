function [pts, weights] = sampleSphere(method, params)

% validate inputs
validateattributes(method, {'char','string'},{},...
    'math.sampleSphere', 'method');

switch method
    case 'lebedev'
        % Lebedev Quadrature method
        % validate inputs
        validateattributes(params, {'numeric'}, {'scalar'},...
                'math.sampleSphere', 'params');
        npts = params(1);
        % Allowed values: {   6,  14,  26,  38,  50,  74,  86, 110, 146,
        %                   170, 194, 230, 266, 302, 350, 434, 590, 770,
        %                   974,1202,1454,1730,2030,2354,2702,3074,3470,
        %                   3890,4334,4802,5294,5810}
        
        samples = math.getLebedevSphere(npts);
        pts = [samples.x, samples.y, samples.z];
        weights = samples.w/(4*pi);
        
    case 'grid'
        % uniform grid on the theta and phi plane
        % validate inputs
        validateattributes(params, {'numeric'}, {'row','ncols',2},...
                'math.sampleSphere', 'params');
        
        nBinsTheta = params(1);
        nBinsPhi = params(2);
        thetaVec = linspace(0, pi, nBinsTheta+1);
        phiVec = linspace(0, 2*pi, nBinsPhi+1);
            
        [thisTheta, thisPhi] = meshgrid(...
        (thetaVec(1:nBinsTheta)+thetaVec(2:nBinsTheta+1))/2,...
        (phiVec(1:nBinsPhi)+phiVec(2:nBinsPhi+1))/2);
        
        x = sin(thisTheta).*cos(thisPhi);
        y = sin(thisTheta).*sin(thisPhi);
        z = cos(thisTheta);
        pts = cat(3, x, y, z);
        
        dS = (thetaVec(2:nBinsTheta+1)-thetaVec(1:nBinsTheta))'*...
             (phiVec(2:nBinsPhi+1)-phiVec(1:nBinsPhi));
        weights = abs(dS/(4*pi)).*sin(thisTheta);
    otherwise
        error('math:sampleSphere',...
        "Unknown method %s. Available methods are 'lebedev' and 'grid'",...
        method);
end
end
            
            