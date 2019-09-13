% Objective: Test if expected value of a numeric function g(r) where r is 
% drawn from a Gamma distribution with shape parameter k and scale v can be
% uniquley attributed to the parametrs k and v; can we estimated paramers k
% and v from E(g(r))?!
%
clc; clear all; close all;

load('/Users/medmfarb/Documents/mycodes/diffusionMRI/caminoInMatlab/test/testData/experimental/scheme.mat');
nScheme = size(scheme, 1);
nBins = 100;

r = linspace(0, 30, nBins+1)*1e-6; % radii in meter
radii = (r(1:nBins) + r(2:nBins+1))*0.5;

% generate signal for each r
g = zeros(nBins, nScheme);
model = cylinder('s0', 1, 'diff', 1e-9, 'theta', pi/3, 'phi', pi/8);
for i = 1:nBins
    model.set('r', radii(i));
    g(i,:) = model.synth(scheme);
end


% compute avarage using gamma distribution
kappa1 = 5;
nu1 = 2.175e-6;
f1 = pdf('gamma', radii, kappa1, nu1);
p1 = diff(cdf('gamma', r, kappa1, nu1));
% scale p for norm 1
p1 = p1./sum(p1);
sig1 = sum(p1'.*g);

% second set of paramters
kappa2 = 10;
nu2 = 1e-6;
f2 = pdf('gamma', radii, kappa2, nu2);
p2 = diff(cdf('gamma', r, kappa2, nu2));
% scale p for norm 1
p2 = p2./sum(p2);
sig2 = sum(p2'.*g);

% ploting
figure; h1 = plot(radii*1e6, f1, 'b--', 'linewidth', 3);
hold on; h2 = plot(radii*1e6, f2, 'r', 'linewidth', 3);
legend([h1, h2], {'k=5, v=2.175e-6', 'k=10, v=1e-6'});
xlabel('r [\mu m]');
title('Gamma Distribution');
set(gca, 'fontsize', 21);

fprintf('Average of signal 1 is %1.6d \n', mean(sig1));
fprintf('Average of signal 2 is %1.6d \n', mean(sig2));
fprintf('RMSE is %1.6d \n', norm(sig1-sig2)/sqrt(nScheme));

figure; scatter(0.5*(sig1+sig2), sig1-sig2, ones(nScheme, 1)*60, 'fill');
title('Bland-Altman Plot');
xlabel('(s1+s2)/2');
ylabel('s1-s2');
set(gca, 'fontsize', 21);
