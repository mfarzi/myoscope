function [long_vectors, radial_vectors, circ_vectors]  = computeLocalRefAxes2(mask, gLongDirection)
    
    if ~iscolumn(gLongDirection)
        gLongDirection = gLongDirection';
    end
    maskSize = size(mask);
    
    idx = find(mask);
    [row, col, dep] = ind2sub(maskSize, idx);
    
    p1 = [mean(col), mean(row), mean(dep)]'; % point to centre
    gLongDirection = gLongDirection/norm(gLongDirection);
    

    p2 = [col, row, dep]; % point to the pixel
    p3 = p1' + ((p2 - p1')*gLongDirection).*gLongDirection';

    % radial vector
    r = p2 - p3;
    r = r ./ sqrt(sum(r.^2, 2));

    % circumferential vector
    c = cross(repmat(gLongDirection', length(r), 1), r);
    c = c ./ sqrt(sum(c.^2, 2));
    
    %% create vector fields
    iLong = zeros(maskSize);
    iLong(mask)= gLongDirection(1);
    jLong = zeros(maskSize);
    jLong(mask)= gLongDirection(2);
    kLong = zeros(maskSize);
    kLong(mask)= gLongDirection(3);
    long_vectors = cat(4, iLong, jLong, kLong);
    
    iRadial = zeros(maskSize);
    iRadial(mask)= r(:,1);
    jRadial = zeros(maskSize);
    jRadial(mask)= r(:,2);
    kRadial = zeros(maskSize);
    kRadial(mask)= r(:,3);
    radial_vectors = cat(4, iRadial, jRadial, kRadial);

    iCirc = zeros(maskSize);
    iCirc(mask)= c(:,1);
    jCirc = zeros(maskSize);
    jCirc(mask)= c(:,2);
    kCirc = zeros(maskSize);
    kCirc(mask)= c(:,3);
    circ_vectors = cat(4, iCirc, jCirc, kCirc);
    
%% test    
thisSlice = false(size(mask));
thisSlice(:, :,16) = mask(:,:,16);
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
%scatter3(y, x, z,'.'); 
axis equal;
hold on; 
quiver3(y, x, z, r2, r1, r3, 0.5, 'r');
%quiver3(y, x, z, l2, l1, l3, 1, 'g');
quiver3(y, x, z, c2, c1, c3, 0.5, 'k');