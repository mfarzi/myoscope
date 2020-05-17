function sheetletAzimuthAngle = computeSheetletAzimuth(ternaryVectorField, circVectors, radialVectors)
%COMPUTEHELIXANGLE Compute Helix angle
%   HA = COMPUTEHELIXANGLE(MASK,)
% Input:
%             vectorField:    Matrix of size 3xN where each column is a
%                             vector
%             longVectors:    Matrix of size 3xN where each column is the
%                             longitudinal vector at one voxel
%             circVectors:    Matrix of size 3xN where each column is the
%                             circumferential vector at one vexel
%
% Email: m.farzi@leeds.ac.uk

% number of voxels
nVoxels = size(ternaryVectorField, 2);
sheetletAzimuthAngle = zeros(nVoxels, 1);

for thisVxl = 1:nVoxels
    % primary eigen vector
    f = ternaryVectorField(:, thisVxl); 
    f = f/norm(f);
    c = circVectors(:, thisVxl);
    c = c/norm(c);
    r = radialVectors(:, thisVxl);
    r = r/norm(r);
    
    % project vectro f onto the cl-plane
    f_cr = dot(c, f) * c + dot(r, f) * r;
    f_cr = f_cr / norm(f_cr);
    
    % compute helix angle
    thisAzimuthAngle = atan2(dot(f_cr, c), dot(f_cr, r));
    
    % since f or -f are the same, helix angle is always between -pi/2 and
    % pi/2
    if thisAzimuthAngle>pi/2
        thisAzimuthAngle = thisAzimuthAngle-pi;
    end

    if thisAzimuthAngle<-pi/2
        thisAzimuthAngle = thisAzimuthAngle+pi;
    end
    
    % convert from radian to degree
    sheetletAzimuthAngle(thisVxl) = thisAzimuthAngle*180/pi;
end