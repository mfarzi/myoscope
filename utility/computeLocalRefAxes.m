function [long_vectors, radial_vectors, circ_vectors]  = computeLocalRefAxes(mask, cavity, gLongVec)
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
% % S1 = abs(X.^2 + Y.^2 - 64 )<9 & (Z<20) & (Z>-20);
% % S2 = abs(X.^2 + Y.^2 - 400 )<30 & (Z<20) & (Z>-20);
% % mask = (X.^2 + Y.^2) > 64 & (X.^2 + Y.^2) < 400 & (Z<20) & (Z>-20);
% % cavity = (X.^2 + Y.^2) <64 & (Z<20) & (Z>-20);
% 
% % sphere
% % S1 = abs(X.^2 + Y.^2 + Z.^2 - 64 )<9;
% % S2 = abs(X.^2 + Y.^2 + Z.^2 - 400 )<30;
% % mask = (X.^2 + Y.^2 + Z.^2) > 64 & (X.^2 + Y.^2 + Z.^2) < 400;
% % cavity = (X.^2 + Y.^2 + Z.^2) <64;
% % volumeViewer(S1|S2);

%% Find the global longitudinal axis by fitting an ellipsoid to points
gLongVec = gLongVec/norm(gLongVec);
     
%%
inner_region = cavity;
outer_region = not(imfill(mask|cavity, 'holes'));

%% solve laplace equcation
% Potential = zeros(size(mask));
% Potential(inner_region) = 0;
% Potential(outer_region) = 10000;
% hFilter = cat(3, [0 0 0; 0 1 0; 0 0 0], [0 1 0; 1 0 1;0 1 0], [0 0 0;0 1 0;0 0 0])/6; 
% hy = cat(3, zeros(3), [0 1 0; 0 0 0;0 -1 0], zeros(3))/2;
% hx = cat(3, zeros(3), [0 1 0; 0 0 0;0 -1 0]', zeros(3))/2;
% hz = cat(3, [0 0 0; 0 1 0; 0 0 0], [0 0 0;0 0 0;0 0 0], [0 0 0;0 -1 0;0 0 0])/2;
% energy_old = 0;
% for iter = 1:500
%     Potential = imfilter(Potential, hFilter);
%     Potential(inner_region) = 0;
%     Potential(outer_region) = 10000;
%     Potential(~mask) = 0;
% 
%     psi_x = imfilter(Potential, hx);
%     psi_y = imfilter(Potential, hy);
%     psi_z = imfilter(Potential, hz);
%     energy_new = sum(sqrt(psi_x(mask).^2+psi_y(mask).^2+psi_z(mask).^2));
%     
%     fprintf('iter %d energy %f\n', iter, (energy_new-energy_old)/energy_old);
%     energy_old = energy_new;
% end

%% The following section comes directly from Jones 2000
% Three-Dimensional Mapping of Cortical Thickness Using Laplace's Equation

% define the potential at each voxel (analogous to voltage)
Potential = zeros(size(mask));
Potential(outer_region) = 10000;
% set the outside to 10000 (arbitrary) and the inside to 0

% initialise the potential field using a distance transform
D_out = bwdist(outer_region);
D_in = bwdist(inner_region);
D = D_out ./ (D_out + D_in);

Potential(mask) = D(mask)* 10000;


% define the structuring element
h = cat(3, [0 0 0; 0 1 0; 0 0 0], [0 1 0; 1 0 1;0 1 0], [0 0 0;0 1 0;0 0 0])/6; %ones([3 3 3])/27;

% run the laplacian until convergence
for its = 1:500 % this should be about 200 according to the paper
    P_plus_one = imfilter(Potential, h, 'replicate');
    
    % only update in the masked region 
    P_plus_one(inner_region) = 0;
    P_plus_one(outer_region) = 10000;
    Potential = P_plus_one;    

end

%% get the direction of the gradient at each voxel
[grad_x, grad_y, grad_z] = gradient(Potential);
gradNorm = sqrt(grad_x.^2 + grad_y.^2 + grad_z.^2);
grad_x = grad_x./gradNorm;
grad_y = grad_y./gradNorm;
grad_z = grad_z./gradNorm;

% store radial vectors
radial_vectors = zeros([size(mask), 3]);
radial_vectors(:,:,:,1) = grad_x;
radial_vectors(:,:,:,2) = grad_y;
radial_vectors(:,:,:,3) = grad_z;
radial_vectors(isnan(radial_vectors)) = 0;
%% 
Global_long = zeros(size(radial_vectors));
Global_long(:,:,:,1) = gLongVec(1);
Global_long(:,:,:,2) = gLongVec(2);
Global_long(:,:,:,3) = gLongVec(3);


% generate the circumferential direction - vector cross product of the
% global longitudinal vectors and the radial vectors
circ_vectors = cross(Global_long, radial_vectors, 4);

% normalise the circumferential vectors (in theory they should already be
% norm 1, but practically I found this wasn't the case)
circ_norm = sqrt(circ_vectors(:,:,:,1).^2 + ...
                 circ_vectors(:,:,:,2).^2 + ...
                 circ_vectors(:,:,:,3).^2);
circ_vectors = circ_vectors ./ repmat(circ_norm, [1 1 1 3]);
circ_vectors(isnan(circ_vectors)) = 0;

% generate the local longitudinal direction
long_vectors = cross(radial_vectors, circ_vectors, 4);

% as before, normalise the vectors
long_norm = sqrt(long_vectors(:,:,:,1).^2 + ...
                 long_vectors(:,:,:,2).^2 + ...
                 long_vectors(:,:,:,3).^2);
long_vectors = long_vectors ./ repmat(long_norm, [1 1 1 3]);
long_vectors(isnan(long_vectors)) = 0;

%% test
thisSlice = false(size(mask));
thisSlice(:,:,17) = mask(:,:,17);
idx = find(thisSlice);
[y, x, z] = ind2sub(size(mask), idx);

clear r1 r2 r3 l1 l2 l3 c1 c2 c3
for i=1:length(x)
    r1(i,1) = radial_vectors(y(i), x(i), z(i), 1);
    r2(i,1) = radial_vectors(y(i), x(i), z(i), 2);
    r3(i,1) = radial_vectors(y(i), x(i), z(i), 3);
end

for i=1:length(x)
    l1(i,1) = long_vectors(y(i), x(i), z(i), 1);
    l2(i,1) = long_vectors(y(i), x(i), z(i), 2);
    l3(i,1) = long_vectors(y(i), x(i), z(i), 3);
end

for i=1:length(x)
    c1(i,1) = circ_vectors(y(i), x(i), z(i), 1);
    c2(i,1) = circ_vectors(y(i), x(i), z(i), 2);
    c3(i,1) = circ_vectors(y(i), x(i), z(i), 3);
end

figure;
scatter3(y, x, z,'.'); 
axis equal;
hold on; 
quiver3(y, x, z, r2, r1, r3, 0, 'r');
%quiver3(y, x, z, l2, l1, l3, 1, 'g');
quiver3(y, x, z, c2, c1, c3, 1, 'k');