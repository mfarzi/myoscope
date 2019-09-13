clc; clear all; close all
addpath('../models');
addpath('../');

load('/Users/medmfarb/Documents/mycodes/diffusionMRI/caminoInMatlab/test/testData/synthetic/scheme.mat');
M = size(scheme, 1);

vFractions1 = [0.42, 0.36, 0.22];
p1a = [1, 5.8e-10, 1.5e-5, 7.6e-6, 0.2, 5.3, 2.3];
p1b = [1, 8.2e-10];
p1c = [1, 9.3e-10, 1.5, 0.8];
model1 = threeCompartmentModel(ellipticalCylinder(p1a), ball(p1b), stick(p1c), 'vFractions', vFractions1);
s1 = model1.synthesize(scheme);

vFractions2 = [0.75, 0, 0.25];
p2a = [1, 7.0e-10, 2.5e-5, 1.4e-5, 0.2, 5.2, 2.4];
p2b = [1, 4.2e-07];
p2c = [1, 8.4e-10, 1.5, 0.8];
model2 = threeCompartmentModel(ellipticalCylinder(p2a), ball(p2b), stick(p2c), 'vFractions', vFractions2);
s2 = model2.synthesize(scheme);

rmse2 = norm(s1-s2)/sqrt(M);

figure('position', [400, 400, 800, 350]); scatter((s2(:)+s1(:))*0.5, s2(:)-s1(:), ones(305,1)*60, 'filled')
set(gca, 'fontsize', 24);
title('Bland-Altman Plot');


model3 = threeCompartmentModel(ellipticalCylinder(), ball(), stick());
model3.fitMultiRun(scheme, s1);
s3 = model3.synthesize(scheme);

rmse3 = norm(s1-s3)/sqrt(M);

figure; scatter((s3(:)+s1(:))*0.5, s3(:)-s1(:));

model4 = ellipticalCylinder();
model4.fitMultiRun(scheme, s1);
s4 = model4.synthesize(scheme);

rmse4 = norm(s1-s4)/sqrt(M);

figure; scatter((s4(:)+s1(:))*0.5, s4(:)-s1(:))

