%clc; close all; clear all

% mu = [0; 0; 1];
% a = 1/2;
% b = 3/2;
% figure;
% for kappa = [4, 8, 16, 32, 64, 128]
% %kappa = 64;
% 
%     c = hypergeom(a, b, kappa);
% 
%     f = @(x, y, z) exp(kappa*(mu(1)*x + mu(2)*y + mu(3)*z).^2);
% 
%     theta = linspace(0, pi*2, 500);
%     phi = pi/2;
%     x = sin(theta)*cos(phi);
%     y = sin(theta)*sin(phi);
%     z = cos(theta);
% 
%     prob = f(x, y, z);
%     prob = prob/max(prob)*250;
% 
%     hold on; plot(y.*prob, z.*prob, 'linewidth', 2);
% end
% 
% % axis equal;
% % xlim([-60, 60]);
% % ylim([-250, 250]);

%%
model = cylinderBD('s0', 1, 'diffPar', 1e-9, 'r', 10e-6, ...
                   'theta', 0, 'phi', 0, 'alpha', 0, ...
                   'nBinsTheta', 200, 'nBinsPhi', 75, ...
                   'kappa1', 4, 'kappa2', 32);

[probs, thisTheta, thisPhi] = model.discritiseBinghamDistribution;

theta = thisTheta(:);
phi = thisPhi(:);
fprintf('integration over unit sphere is: %1.2f\n', sum(probs.*sin(theta)));


x = sin(theta).*cos(phi);
y = sin(theta).*sin(phi);
z = cos(theta);

P = reshape(probs, size(thisTheta));
X = sin(thisTheta).*cos(thisPhi).*P;
Y = sin(thisTheta).*sin(thisPhi).*P;
Z = cos(thisTheta).*P;

scale = max(Z(:));
X = X/scale*8;
Y = Y/scale*8;
Z = Z/scale*8;


figure('position', [400 400 400 900]);
mesh(X, Y, Z, P);
colormap('jet');
daspect([1 1 1])
view(-50,25);
xlim([-2 2]);
ylim([-2 2]);
zlim([-8 8]);