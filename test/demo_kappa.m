% Objective: Test if expected value of a numeric function g(r) where r is 
% drawn from a Gamma distribution with shape parameter k and scale v can be
% uniquley attributed to the parametrs k and v; can we estimated paramers k
% and v from E(g(r))?!
%
clc; clear all; close all;

nBins = 100;
r = linspace(0, 30, nBins+1)*1e-6; % radii in meter
radii = (r(1:nBins) + r(2:nBins+1))*0.5;

% compute avarage using gamma distribution
nu = 1.5e-6;
figure;
hold on;
for i = 1:10
    kappa = i;
    f = pdf('gamma', radii, kappa, nu);
    plot(radii*1e6, f, 'linewidth', 3);
end

xlabel('r [\mu m]');
title('Gamma Distribution');
set(gca, 'fontsize', 21);
