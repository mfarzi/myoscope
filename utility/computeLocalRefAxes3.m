function [long_vectors, radial_vectors, circ_vectors]  = computeLocalRefAxes3(mask, Sbase, Sapex, Sepi, Slv, Srv)
% COMPUTELOCALREFAXES Extract local longitudinal, circumferential, and 
% radial axes from a binary mask
%
% % Email:m.farzi@leeds.ac.uk

%% extract two surfaces
% nRows = 100;
% nCols = 100;
% nDeps = 100;
% [X, Y, Z] = meshgrid(-(nCols-1)/2:(nCols-1)/2, ...
%                      -(nRows-1)/2:(nRows-1)/2, ...
%                      -(nDeps-1)/2:(nDeps-1)/2);
% 
% % cylinder                 
% S1 = abs(X.^2 + Y.^2 - 64 )<9 & (Z<20) & (Z>-20);
% S2 = abs(X.^2 + Y.^2 - 400 )<30 & (Z<20) & (Z>-20);
% S3 = (X.^2 + Y.^2) > 64 & (X.^2 + Y.^2) < 400 & abs(Z-19.5)<0.5;
% S3p = (X.^2 + Y.^2) > 64 & (X.^2 + Y.^2) < 400 & abs(Z-18.5)<0.5;
% S4 = (X.^2 + Y.^2) > 64 & (X.^2 + Y.^2) < 400 & abs(Z+19.5)<0.5;
% S4p = (X.^2 + Y.^2) > 64 & (X.^2 + Y.^2) < 400 & abs(Z+18.5)<0.5;

%[f, v] = isosurface(X, Y, Z, S3);
%n = isonormals(X,Y,Z,S3,v);
%figure; quiver3(v(:,1), v(:,2), v(:,3), n(:,1), n(:,2), n(:,3));

%mask = (X.^2 + Y.^2) > 64 & (X.^2 + Y.^2) < 400 & (Z<20) & (Z>-20);
%cavity = (X.^2 + Y.^2) <64 & (Z<20) & (Z>-20);

% sphere
% S1 = abs(X.^2 + Y.^2 + Z.^2 - 64 )<9;
% S2 = abs(X.^2 + Y.^2 + Z.^2 - 400 )<30;
% mask = (X.^2 + Y.^2 + Z.^2) > 64 & (X.^2 + Y.^2 + Z.^2) < 400;
% cavity = (X.^2 + Y.^2 + Z.^2) <64;
%volumeViewer(S1|S4|S4p);

%% compute 
potential_radial = solveLaplace(mask, Sepi, 0, Slv, 100, Srv, 10, 500);
potential_long   = solveLaplace_long(mask, Sbase, 0, Sapex, 100, 500);     
%% get the direction of the gradient at each voxel
r = grad(potential_radial, mask);
l = grad(potential_long, mask);

rdotl = sum(r.*l, 2);
r = r - rdotl.*l;
r = r./sqrt(sum(r.^2,2));

%c = cross(l,r, 2);

% store radial vectors
radial_vectors_x = zeros(size(mask));
radial_vectors_x(mask) = r(:,1);
radial_vectors_y = zeros(size(mask));
radial_vectors_y(mask) = r(:,2);
radial_vectors_z = zeros(size(mask));
radial_vectors_z(mask) = r(:,3);
radial_vectors = cat(4, radial_vectors_x, radial_vectors_y, radial_vectors_z);
radial_vectors(isnan(radial_vectors)) = 0;

% store radial vectors
long_vectors_x = zeros(size(mask));
long_vectors_x(mask) = l(:,1);
long_vectors_y = zeros(size(mask));
long_vectors_y(mask) = l(:,2);
long_vectors_z = zeros(size(mask));
long_vectors_z(mask) = l(:,3);
long_vectors = cat(4, long_vectors_x, long_vectors_y, long_vectors_z);
long_vectors(isnan(long_vectors)) = 0;

circ_vectors = cross(long_vectors, radial_vectors, 4);
circ_vectors(isnan(circ_vectors)) = 0;

%% test
% thisSlice = false(size(mask));
% thisSlice(:,:,17) = mask(:,:,17);
% idx = find(thisSlice);
% [y, x, z] = ind2sub(size(mask), idx);
% 
% clear r1 r2 r3 l1 l2 l3 c1 c2 c3
% for i=1:length(x)
%     r1(i,1) = radial_vectors(y(i), x(i), z(i), 1);
%     r2(i,1) = radial_vectors(y(i), x(i), z(i), 2);
%     r3(i,1) = radial_vectors(y(i), x(i), z(i), 3);
% end
% 
% for i=1:length(x)
%     l1(i,1) = long_vectors(y(i), x(i), z(i), 1);
%     l2(i,1) = long_vectors(y(i), x(i), z(i), 2);
%     l3(i,1) = long_vectors(y(i), x(i), z(i), 3);
% end
% 
% for i=1:length(x)
%     c1(i,1) = circ_vectors(y(i), x(i), z(i), 1);
%     c2(i,1) = circ_vectors(y(i), x(i), z(i), 2);
%     c3(i,1) = circ_vectors(y(i), x(i), z(i), 3);
% end
% 
% figure;
% scatter3(y, x, z,'.'); 
% axis equal;
% hold on; 
% quiver3(y, x, z, r2, r1, r3, 0, 'r');
% quiver3(y, x, z, l2, l1, l3, 1, 'g');
% quiver3(y, x, z, c2, c1, c3, 1, 'k');
end

function Potential = solveLaplace(mask, S1, phi1, S2, phi2, S3, phi3, maxIter)
    idx = find(mask);
    nPxl = length(idx);
    [row, col, dep] = ind2sub(size(mask), idx);

    Potential = zeros(size(mask));
    Potential(S1) = phi1;
    Potential(S2) = phi2;
    Potential(S3) = phi3;
    for iter = 1:maxIter
        p = arrayfun(@(i) updatePotential(Potential, mask, row(i), col(i), dep(i)), 1:nPxl);
        Potential(mask) = p;
        Potential(~mask) = 0;
        Potential(S1) = phi1;
        Potential(S2) = phi2;
        Potential(S3) = phi3;
    end
end

function Potential = solveLaplace_long(mask, S1, phi1, S2, phi2, maxIter)
    idx = find(mask);
    nPxl = length(idx);
    [row, col, dep] = ind2sub(size(mask), idx);

    Potential = zeros(size(mask));
    Potential(S1) = phi1;
    Potential(S2) = phi2;
    for iter = 1:maxIter
        p = arrayfun(@(i) updatePotential(Potential, mask, row(i), col(i), dep(i)), 1:nPxl);
        Potential(mask) = p;
        Potential(~mask) = 0;
        Potential(S1) = phi1;
        Potential(S2) = phi2;
    end
end

function val = updatePotential(Potential, mask, row, col, dep)
    xflag = false;
    yflag = false;
    zflag = false;
    
    if mask(row+1, col, dep) && mask(row-1, col, dep)
        y = Potential(row+1, col, dep) + Potential(row-1, col, dep);
    elseif mask(row+1, col, dep) && ~mask(row-1, col, dep)
        y = Potential(row+1, col, dep)*2;
    elseif ~mask(row+1, col, dep) && mask(row-1, col, dep)
        y = Potential(row-1, col, dep)*2;
    else
        y = 0;
        yflag = true;
    end
    
    if mask(row, col+1, dep) && mask(row, col-1, dep)
        x = Potential(row, col+1, dep) + Potential(row, col-1, dep);
    elseif mask(row, col+1, dep) && ~mask(row, col-1, dep)
        x = Potential(row, col+1, dep)*2;
    elseif ~mask(row, col+1, dep) && mask(row, col-1, dep)
        x = Potential(row, col-1, dep)*2;
    else
        x = 0;
        xflag = true;
    end
    
    if mask(row, col, dep+1) && mask(row, col, dep-1)
        z = Potential(row, col, dep+1) + Potential(row, col, dep-1);
    elseif mask(row, col, dep+1) && ~mask(row, col, dep-1)
        z = Potential(row, col, dep+1)*2;
    elseif ~mask(row, col, dep+1) && mask(row, col, dep-1)
        z = Potential(row, col, dep-1)*2;
    else
        z = 0;
        zflag = true;
    end
    
    val = (x+y+z)/6;
    
    if xflag&&yflag&&zflag
        fprintf('[%d %d %d]->%d\n', row, col, dep, 1);
    end
        
    
end

function gx = gradx(Potential, mask, row, col, dep)
    if mask(row, col+1, dep) && mask(row, col-1, dep)
        gx = (Potential(row, col-1, dep) - Potential(row, col+1, dep))/2;
    elseif mask(row, col+1, dep) && ~mask(row, col-1, dep)
        gx = Potential(row, col, dep)-Potential(row, col+1, dep);
    elseif ~mask(row, col+1, dep) && mask(row, col-1, dep)
        gx = Potential(row, col-1, dep) - Potential(row, col, dep);
    else
        gx = 1;
    end
end

function gy = grady(Potential, mask, row, col, dep)
    if mask(row+1, col, dep) && mask(row-1, col, dep)
        gy = (Potential(row-1, col, dep) -Potential(row+1, col, dep))/2;
    elseif mask(row+1, col, dep) && ~mask(row-1, col, dep)
        gy = Potential(row, col, dep)-Potential(row+1, col, dep);
    elseif ~mask(row+1, col, dep) && mask(row-1, col, dep)
        gy = Potential(row-1, col, dep) - Potential(row, col, dep);
    else
        gy = 1;
    end
end

function gz = gradz(Potential, mask, row, col, dep)
    if mask(row, col, dep+1) && mask(row, col, dep-1)
        gz = (Potential(row, col, dep-1) -Potential(row, col, dep+1))/2;
    elseif mask(row, col, dep+1) && ~mask(row, col, dep-1)
        gz = Potential(row, col, dep)-Potential(row, col, dep+1);
    elseif ~mask(row, col, dep+1) && mask(row, col, dep-1)
        gz = Potential(row, col, dep-1) - Potential(row, col, dep);
    else
        gz = 1;
    end
end

function g = grad(potential_radial, mask)
    idx = find(mask);
    nPxl = length(idx);
    [row, col, dep] = ind2sub(size(mask), idx);
    
    grad_x = arrayfun(@(i) gradx(potential_radial, mask, row(i), col(i), dep(i)), 1:nPxl);
    grad_y = arrayfun(@(i) grady(potential_radial, mask, row(i), col(i), dep(i)), 1:nPxl);
    grad_z = arrayfun(@(i) gradz(potential_radial, mask, row(i), col(i), dep(i)), 1:nPxl);
    gradNorm = sqrt(grad_x.^2 + grad_y.^2 + grad_z.^2);
    grad_x = grad_x./gradNorm;
    grad_y = grad_y./gradNorm;
    grad_z = grad_z./gradNorm;
    g = [grad_x', grad_y', grad_z'];
end
