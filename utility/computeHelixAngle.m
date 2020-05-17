function helixAngle = computeHelixAngle(vectorField, longVectors, circVectors)
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
nVoxels = size(vectorField, 2);
helixAngle = zeros(nVoxels, 1);

for thisVxl = 1:nVoxels
    % primary eigen vector
    f = vectorField(:, thisVxl); 
    f = f/norm(f);
    c = circVectors(:, thisVxl);
    c = c/norm(c);
    l = longVectors(:, thisVxl);
    l = l/norm(l);
    
    % project vectro f onto the cl-plane
    f_cl = dot(c, f) * c + dot(l, f) * l;
    f_cl = f_cl / norm(f_cl);
    
    % compute helix angle
    thisHelixAngle = atan2(dot(f_cl, l), dot(f_cl, c));
    
    % since f or -f are the same, helix angle is always between -pi/2 and
    % pi/2
    if thisHelixAngle>pi/2
        thisHelixAngle = thisHelixAngle-pi;
    end

    if thisHelixAngle<-pi/2
        thisHelixAngle = thisHelixAngle+pi;
    end
    
    % convert from radian to degree
    helixAngle(thisVxl) = thisHelixAngle*180/pi;
end