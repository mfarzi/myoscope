function sheetletElevationAngle = computeSheetletElevation(ternaryVectorField, longVectors, radialVectors)
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
sheetletElevationAngle = zeros(nVoxels, 1);

for thisVxl = 1:nVoxels
    % primary eigen vector
    f = ternaryVectorField(:, thisVxl); 
    f = f/norm(f);
    l = longVectors(:, thisVxl);
    l = l/norm(l);
    r = radialVectors(:, thisVxl);
    r = r/norm(r);
    
    % project vectro f onto the cl-plane
    f_rl = dot(l, f) * l + dot(r, f) * r;
    f_rl = f_rl / norm(f_rl);
    
    % compute helix angle
    thisElevationAngle = atan2(dot(f_rl, l), dot(f_rl, r));
    
    % since f or -f are the same, helix angle is always between -pi/2 and
    % pi/2
    if thisElevationAngle>pi/2
        thisElevationAngle = thisElevationAngle-pi;
    end

    if thisElevationAngle<-pi/2
        thisElevationAngle = thisElevationAngle+pi;
    end
    
    % convert from radian to degree
    sheetletElevationAngle(thisVxl) = thisElevationAngle*180/pi;
end