% This script test if each model synthesizes the signal correctly or not
clc; clear all; close all
addpath('../models');
addpath('../linker');
addpath('../'); 

path_to_data = ['/Users/medmfarb/Documents/mycodes/diffusionMRI/',...
                'caminoInMatlab/test/testData/synthetic/'];
            
tol = 1e-6;                     % This error is due to machine precision
                                % for camino (java version)            
%% tensor
load(fullfile(path_to_data, 'tensor.mat'), 'testData', 'scheme');
N = length(testData);
nScheme = size(scheme, 1);

rmse = zeros(N, 1);
s1 = [];
s2 = [];
for n = 1:N
    params = [testData(n).diffPar, testData(n).diffPerp1, testData(n).diffPerp2, ...
              testData(n).theta, testData(n).phi, testData(n).alpha];
    model = tensor([1, params]);
    s = model.synthesize(scheme);
    rmse(n) = norm(s-testData(n).s)/sqrt(nScheme);
    
    s1 = [s1; s];
    s2 = [s2; testData(n).s];
end

if max(rmse)>tol
    fprintf('Validation test is faild for class %s (max rmse = %1.2e).\n', model.name, max(rmse));
else
    fprintf('Validation test is passed for class %s (max rmse = %1.2e).\n', model.name, max(rmse));
end

figure('position', [400 400 900 350]); scatter((s1+s2)/2, s1-s2, ones(length(s1),1)*50, 'filled');
title({'Bland-Altman Plot', sprintf('( %s )', model.name)}, 'fontsize', 24);
set(gca, 'fontsize', 22);

clear test scheme model rms_err

%% ball
load(fullfile(path_to_data, 'ball.mat'), 'testData', 'scheme');
N = length(testData);
nScheme = size(scheme, 1);

rmse = zeros(N, 1);
s1 = [];
s2 = [];
for n = 1:N
    model = ball([1; testData(n).diff]);
    s = model.synthesize(scheme);
    rmse(n) = norm(s-testData(n).s)/sqrt(nScheme);
    
    s1 = [s1; s];
    s2 = [s2; testData(n).s];
end

if max(rmse)>tol
    fprintf('Validation test is faild for class %s (max rmse = %1.2e).\n', model.name, max(rmse));
else
    fprintf('Validation test is passed for class %s (max rmse = %1.2e).\n', model.name, max(rmse));
end

figure('position', [400 400 900 350]); scatter((s1+s2)/2, s1-s2, ones(length(s1),1)*50, 'filled');
title({'Bland-Altman Plot', sprintf('( %s )', model.name)}, 'fontsize', 24);
set(gca, 'fontsize', 22);

clear testData scheme model rms_err

%% zeppelin
load(fullfile(path_to_data, 'zeppelin.mat'), 'testData', 'scheme');
N = length(testData);
nScheme = size(scheme, 1);

rmse = zeros(N, 1);
s1 = [];
s2 = [];
for n = 1:N
    model = zeppelin([1, testData(n).diffPar, testData(n).diffPerp, ...
                      testData(n).theta, testData(n).phi]);
    s = model.synthesize(scheme);
    rmse(n) = norm(s-testData(n).s)/sqrt(nScheme);
    
    s1 = [s1; s];
    s2 = [s2; testData(n).s];
end

if max(rmse)>tol
    fprintf('Validation test is faild for class %s (max rmse = %1.2e).\n', model.name, max(rmse));
else
    fprintf('Validation test is passed for class %s (max rmse = %1.2e).\n', model.name, max(rmse));
end

figure('position', [400 400 900 350]); scatter((s1+s2)/2, s1-s2, ones(length(s1),1)*50, 'filled');
title({'Bland-Altman Plot', sprintf('( %s )', model.name)}, 'fontsize', 24);
set(gca, 'fontsize', 22);

clear testData scheme model rms_err

%% stick
load(fullfile(path_to_data, 'stick.mat'), 'testData', 'scheme');
N = length(testData);
nScheme = size(scheme, 1);

rmse = zeros(N, 1);
s1 = [];
s2 = [];
for n = 1:N
    s0 = 1;
    params = [s0, testData(n).diff, testData(n).theta, testData(n).phi];
    model = stick(params);
    s = model.synthesize(scheme);
    rmse(n) = norm(s-testData(n).s)/sqrt(nScheme);
    
    s1 = [s1; s];
    s2 = [s2; testData(n).s];
end

if max(rmse)>tol
    fprintf('Validation test is faild for class %s (max rmse = %1.2e).\n', model.name, max(rmse));
else
    fprintf('Validation test is passed for class %s (max rmse = %1.2e).\n', model.name, max(rmse));
end

figure('position', [400 400 900 350]); scatter((s1+s2)/2, s1-s2, ones(length(s1),1)*50, 'filled');
title({'Bland-Altman Plot', sprintf('( %s )', model.name)}, 'fontsize', 24);
set(gca, 'fontsize', 22);

clear testData scheme model rms_err

%% cylinderGPD
load(fullfile(path_to_data, 'cylinderGPD.mat'), 'testData', 'scheme');
N = length(testData);
nScheme = size(scheme, 1);

rmse = zeros(N, 1);
s1 = [];
s2 = [];
for n = 1:N
    s0 = 1;
    params = [s0, testData(n).diff, testData(n).r, testData(n).theta, testData(n).phi];
    model = cylinderGPD(params);
    s = model.synth(scheme);
    rmse(n) = norm(s-testData(n).s)/sqrt(nScheme);
    
    s1 = [s1; s];
    s2 = [s2; testData(n).s];
end

if max(rmse)>tol
    fprintf('Validation test is faild for class %s (max rmse = %1.2e).\n', model.name, max(rmse));
else
    fprintf('Validation test is passed for class %s (max rmse = %1.2e).\n', model.name, max(rmse));
end

figure('position', [400 400 900 350]); scatter((s1+s2)/2, s1-s2, ones(length(s1),1)*50, 'filled');
title({'Bland-Altman Plot', sprintf('( %s )', model.name)}, 'fontsize', 24);
set(gca, 'fontsize', 22);
clear testData scheme model rms_err

%% cylinderGDR
load(fullfile(path_to_data, 'cylinderGDR.mat'), 'testData', 'scheme');
N = length(testData);
nScheme = size(scheme, 1);

rmse = zeros(N, 1);
s1 = [];
s2 = [];
for n = 1:N
    s0 = 1;
    params = [s0, testData(n).diff, testData(n).kappa, testData(n).nu, testData(n).theta, testData(n).phi];
    model = cylinderGDR(params');
    % NOTE THAT CAMINO-IN-JAVA USE NBINS=5
    model.set('nbins', 5);
    s = model.synth(scheme);
    rmse(n) = norm(s-testData(n).s)/sqrt(nScheme);
    
    s1 = [s1; s];
    s2 = [s2; testData(n).s];
end

if max(rmse)>tol
    fprintf('Validation test is faild for class %s (max rmse = %1.2e).\n', model.name, max(rmse));
else
    fprintf('Validation test is passed for class %s (max rmse = %1.2e).\n', model.name, max(rmse));
end

figure('position', [400 400 900 350]); scatter((s1+s2)/2, s1-s2, ones(length(s1),1)*50, 'filled');
title({'Bland-Altman Plot', sprintf('( %s )', model.name)}, 'fontsize', 24);
set(gca, 'fontsize', 22);
% clear testData scheme model rms_err

%% ellipticalCylinder
load(fullfile(path_to_data, 'ellipticalCylinder.mat'), 'testData', 'scheme');
N = length(testData);
nScheme = size(scheme, 1);

rmse = zeros(N, 1);
s1 = [];
s2 = [];
for n = 1:N
    s0 = 1;
    params = [s0, testData(n).diff, testData(n).r1, testData(n).r2, ...
                  testData(n).theta, testData(n).phi, testData(n).alpha];
    model = ellipticalCylinder(params);
    s = model.synthesize(scheme);
    rmse(n) = norm(s-testData(n).s)/sqrt(nScheme);
    
    s1 = [s1; s];
    s2 = [s2; testData(n).s];
end

if max(rmse)>tol
    fprintf('Validation test is faild for class %s (max rmse = %1.2e).\n', model.name, max(rmse));
else
    fprintf('Validation test is passed for class %s (max rmse = %1.2e).\n', model.name, max(rmse));
end

figure('position', [400 400 900 350]); scatter((s1+s2)/2, s1-s2, ones(length(s1),1)*50, 'filled');
title({'Bland-Altman Plot', sprintf('( %s )', model.name)}, 'fontsize', 24);
set(gca, 'fontsize', 22);

clear testData scheme model rms_err

%% ball_stick
load(fullfile(path_to_data, 'ball_stick.mat'), 'testData', 'scheme');
N = length(testData);
nScheme = size(scheme, 1);

rmse = zeros(N, 1);
s1 = [];
s2 = [];
for n = 1:N
    s0 = 1;
    comp1 = ball([s0, testData(n).comp1_diff]);
    comp2 = stick([s0, testData(n).comp2_diff, testData(n).comp2_theta, testData(n).comp2_phi]);
    vFractions = [testData(n).comp1_vFraction, testData(n).comp2_vFraction];
    model = twoCompartmentModel(comp1, comp2, 'vFractions', vFractions);
    s = model.synthesize(scheme);
    rmse(n) = norm(s-testData(n).s)/sqrt(nScheme);
    
    s1 = [s1; s];
    s2 = [s2; testData(n).s];
end

if max(rmse)>tol
    fprintf('Validation test is faild for class %s (max rmse = %1.2e).\n', model.name, max(rmse));
else
    fprintf('Validation test is passed for class %s (max rmse = %1.2e).\n', model.name, max(rmse));
end

figure('position', [400 400 900 350]); scatter((s1+s2)/2, s1-s2, ones(length(s1),1)*50, 'filled');
title({'Bland-Altman Plot', sprintf('( %s-%s )', comp1.name, comp2.name)}, 'fontsize', 24);
set(gca, 'fontsize', 22);

clear testData scheme model rms_err

%% cyllinder_zeppelin
load(fullfile(path_to_data, 'cylinder_zeppelin.mat'), 'testData', 'scheme');
N = length(testData);
nScheme = size(scheme, 1);

rmse = zeros(N, 1);
s1 = [];
s2 = [];
for n = 1:N
    s0 = 1;
    comp1 = cylinder([s0, testData(n).comp1_diff, testData(n).comp1_r, ...
                      testData(n).comp1_theta, testData(n).comp1_phi]);
    comp2 = zeppelin([s0, testData(n).comp2_diffPar, ...
                      testData(n).comp2_diffPerp, testData(n).comp2_theta,...
                      testData(n).comp2_phi]);
    vFractions = [testData(n).comp1_vFraction, testData(n).comp2_vFraction];
    model = twoCompartmentModel(comp1, comp2, 'vFractions', vFractions);
    s = model.synthesize(scheme);
    rmse(n) = norm(s-testData(n).s)/sqrt(nScheme);
    
    s1 = [s1; s];
    s2 = [s2; testData(n).s];
end

if max(rmse)>tol
    fprintf('Validation test is faild for class %s (max rmse = %1.2e).\n', model.name, max(rmse));
else
    fprintf('Validation test is passed for class %s (max rmse = %1.2e).\n', model.name, max(rmse));
end

figure('position', [400 400 900 350]); scatter((s1+s2)/2, s1-s2, ones(length(s1),1)*50, 'filled');
title({'Bland-Altman Plot', sprintf('( %s-%s )', comp1.name, comp2.name)}, 'fontsize', 24);
set(gca, 'fontsize', 22);

clear testData scheme model rms_err

%% ellitpicalCylinder_ball_stick
load(fullfile(path_to_data, 'ellipticalCylinder_ball_stick.mat'), 'testData', 'scheme');
N = length(testData);
nScheme = size(scheme, 1);

rmse = zeros(N, 1);
s1 = [];
s2 = [];
for n = 1:N
    s0 = 1;
    comp1 = ellipticalCylinder([s0, testData(n).comp1_diff, testData(n).comp1_r1, testData(n).comp1_r2, testData(n).comp1_theta, testData(n).comp1_phi, testData(n).comp1_alpha]);
    comp2 = ball([s0, testData(n).comp2_diff]);
    comp3 = stick([s0, testData(n).comp3_diff, testData(n).comp3_theta, testData(n).comp3_phi]);
    vFractions = [testData(n).comp1_vFraction, testData(n).comp2_vFraction, testData(n).comp3_vFraction];
    model = threeCompartmentModel(comp1, comp2, comp3, 'vFractions', vFractions);
    s = model.synthesize(scheme);
    rmse(n) = norm(s-testData(n).s)/sqrt(nScheme);
    
    s1 = [s1; s];
    s2 = [s2; testData(n).s];
end

if max(rmse)>tol
    fprintf('Validation test is faild for class %s (max rmse = %1.2e).\n', model.name, max(rmse));
else
    fprintf('Validation test is passed for class %s (max rmse = %1.2e).\n', model.name, max(rmse));
end

figure('position', [400 400 900 350]); scatter((s1+s2)/2, s1-s2, ones(length(s1),1)*50, 'filled');
title({'Bland-Altman Plot', sprintf('( %s-%s-%s )', comp1.name, comp2.name, comp3.name)}, 'fontsize', 24);
set(gca, 'fontsize', 22);

clear testData scheme model rms_err