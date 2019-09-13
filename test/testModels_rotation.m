% This script test if each compartent model has implemented the method
% 'rotateAxis' properly to remove ambiguity with theta-phi-alpha
% presentation in camino.
clc; clear all; close all
addpath('../models');
addpath('../');

path_to_data = ['/Users/medmfarb/Documents/mycodes/diffusionMRI/',...
                'caminoInMatlab/test/testData/synthetic/'];

tol = 1e-12;            
%% tensor         
load(fullfile(path_to_data, 'tensor.mat'), 'testData', 'scheme');
N = length(testData);
nScheme = size(scheme, 1);

rmse = zeros(N, 1);
for n = 1:100
    s0 = 1;
    params = [s0, testData(n).diffPar, testData(n).diffPerp1, testData(n).diffPerp2, ...
              testData(n).theta, testData(n).phi, testData(n).alpha];
    model = tensor(params);  
    s1 = model.synthesize(scheme);
    
    % rotate the axis and synthesize signals again
    model.rotateAxis();
    s2 = model.synthesize(scheme);
    
    rmse(n) = norm(s1-s2)/sqrt(nScheme);
    
    % correct parameters
    p = model.getParams('excludeS0', true);
    testData(n,1).diffPar   = p(1);
    testData(n,1).diffPerp1 = p(2);
    testData(n,1).diffPerp2 = p(3);
    testData(n,1).theta     = p(4);
    testData(n,1).phi       = p(5);
    testData(n,1).alpha     = p(6);
end            

if max(rmse) < tol
    fprintf('rotateAxis function works fine in class %s (max rmse = %1.2e).\n', model.name, max(rmse));
    save('/Users/medmfarb/Documents/mycodes/diffusionMRI/caminoInMatlab/test/testData/synthetic/tensor_rotated.mat', 'testData', 'scheme');
else
    error('rotateAxis function is buggy in class %s! (max rmse = %1.2e)\n', model.name, max(rmse));
end

%% ball         
load(fullfile(path_to_data, 'ball.mat'), 'testData', 'scheme');
N = length(testData);
nScheme = size(scheme, 1);

rmse = zeros(N, 1);
for n = 1:100
    s0 = 1;
    params = [s0, testData(n).diff];
    model = ball(params);  
    s1 = model.synthesize(scheme);
    
    % rotate the axis and synthesize signals again
    model.rotateAxis();
    s2 = model.synthesize(scheme);
    
    rmse(n) = norm(s1-s2)/sqrt(nScheme);
    
    % correct parameters
    p = model.getParams('excludeS0', true);
    testData(n,1).diff = p(1);
end            

if max(rmse) < tol
    fprintf('rotateAxis function works fine in class %s (max rmse = %1.2e).\n', model.name, max(rmse));
    save('/Users/medmfarb/Documents/mycodes/diffusionMRI/caminoInMatlab/test/testData/synthetic/ball_rotated.mat', 'testData', 'scheme');
else
    error('rotateAxis function is buggy in class %s! (max rmse = %1.2e)\n', model.name, max(rmse));
end

%% zeppelin         
load(fullfile(path_to_data, 'zeppelin.mat'), 'testData', 'scheme');
N = length(testData);
nScheme = size(scheme, 1);

rmse = zeros(N, 1);
for n = 1:100
    s0 = 1;
    params = [s0, testData(n).diffPar, testData(n).diffPerp, ...
              testData(n).theta, testData(n).phi];
    model = zeppelin(params);  
    s1 = model.synthesize(scheme);
    
    % rotate the axis and synthesize signals again
    model.rotateAxis();
    s2 = model.synthesize(scheme);
    
    rmse(n) = norm(s1-s2)/sqrt(nScheme);
    
    % correct parameters
    p = model.getParams('excludeS0', true);
    testData(n,1).diffPar  = p(1);
    testData(n,1).diffPerp = p(2);
    testData(n,1).theta    = p(3);
    testData(n,1).phi      = p(4);
end            

if max(rmse) < tol
    fprintf('rotateAxis function works fine in class %s (max rmse = %1.2e).\n', model.name, max(rmse));
    save('/Users/medmfarb/Documents/mycodes/diffusionMRI/caminoInMatlab/test/testData/synthetic/zeppelin_rotated.mat', 'testData', 'scheme');
else
    error('rotateAxis function is buggy in class %s! (max rmse = %1.2e)\n', model.name, max(rmse));
end

%% stick       
load(fullfile(path_to_data, 'stick.mat'), 'testData', 'scheme');
N = length(testData);
nScheme = size(scheme, 1);

rmse = zeros(N, 1);
for n = 1:100
    s0 = 1;
    params = [s0, testData(n).diff, testData(n).theta, testData(n).phi];
    model = stick(params);  
    s1 = model.synthesize(scheme);
    
    % rotate the axis and synthesize signals again
    model.rotateAxis();
    s2 = model.synthesize(scheme);
    
    rmse(n) = norm(s1-s2)/sqrt(nScheme);
    
    % correct parameters
    p = model.getParams('excludeS0', true);
    testData(n,1).diff  = p(1);
    testData(n,1).theta = p(2);
    testData(n,1).phi   = p(3);
end            

if max(rmse) < tol
    fprintf('rotateAxis function works fine in class %s (max rmse = %1.2e).\n', model.name, max(rmse));
    save('/Users/medmfarb/Documents/mycodes/diffusionMRI/caminoInMatlab/test/testData/synthetic/stick_rotated.mat', 'testData', 'scheme');
else
    error('rotateAxis function is buggy in class %s! (max rmse = %1.2e)\n', model.name, max(rmse));
end

%% cylinder
load(fullfile(path_to_data, 'cylinder.mat'), 'testData', 'scheme');
N = length(testData);
nScheme = size(scheme, 1);

rmse = zeros(N, 1);
for n = 1:100
    s0 = 1;
    params = [s0, testData(n).diff, testData(n).r, testData(n).theta, ...
              testData(n).phi];
    model = cylinder(params);  
    s1 = model.synthesize(scheme);
    
    % rotate the axis and synthesize signals again
    model.rotateAxis();
    s2 = model.synthesize(scheme);
    
    rmse(n) = norm(s1-s2)/sqrt(nScheme);
    
    % correct parameters
    p = model.getParams('excludeS0', true);
    testData(n,1).diff  = p(1);
    testData(n,1).r     = p(2);
    testData(n,1).theta = p(3);
    testData(n,1).phi   = p(4);
    
end            

if max(rmse) < tol
    fprintf('rotateAxis function works fine in class %s (max rmse = %1.2e).\n', model.name, max(rmse));
    save('/Users/medmfarb/Documents/mycodes/diffusionMRI/caminoInMatlab/test/testData/synthetic/cylinder_rotated.mat', 'testData', 'scheme');
else
    error('rotateAxis function is buggy in class %s! (max rmse = %1.2e)\n', model.name, max(rmse));
end

%% ellipticalyCylinder
load(fullfile(path_to_data, 'ellipticalCylinder.mat'), 'testData', 'scheme');
N = length(testData);
nScheme = size(scheme, 1);

rmse = zeros(N, 1);
for n = 1:100
    s0 = 1;
    params = [s0, testData(n).diff, testData(n).r1, testData(n).r2, ...
              testData(n).theta, testData(n).phi, testData(n).alpha];
    model = ellipticalCylinder(params);  
    s1 = model.synthesize(scheme);
    
    % rotate the axis and synthesize signals again
    model.rotateAxis();
    s2 = model.synthesize(scheme);
    
    rmse(n) = norm(s1-s2)/sqrt(nScheme);
    
    % correct parameters
    p = model.getParams('excludeS0', true);
    testData(n,1).diff  = p(1);
    testData(n,1).r1    = p(2);
    testData(n,1).r2    = p(3);
    testData(n,1).theta = p(4);
    testData(n,1).phi   = p(5);
    testData(n,1).alpha = p(6);
    
end            

if max(rmse) < tol
    fprintf('rotateAxis function works fine in class %s (max rmse = %1.2e).\n', model.name, max(rmse));
    save('/Users/medmfarb/Documents/mycodes/diffusionMRI/caminoInMatlab/test/testData/synthetic/ellipticalCylinder_rotated.mat', 'testData', 'scheme');
else
     error('rotateAxis function is buggy in class %s! (max rmse = %1.2e)\n', model.name, max(rmse));
end

%% ball_stick
load(fullfile(path_to_data, 'ball_stick.mat'), 'testData', 'scheme');
N = length(testData);
nScheme = size(scheme, 1);

rmse = zeros(N, 1);
for n = 1:100
    s0 = 1;
    
    % read the input parameters
    comp1 = ball([s0, testData(n).comp1_diff]);
    comp2 = stick([s0, testData(n).comp2_diff, testData(n).comp2_theta, ...
                   testData(n).comp2_phi]);
    
    vFractions = [testData(n).comp1_vFraction, testData(n).comp2_vFraction];
    model = twoCompartmentModel(comp1, comp2, 'vFractions', vFractions);  
    s1 = model.synthesize(scheme);
    
    % rotate the axis and synthesize signals again
    model.rotateAxis();
    s2 = model.synthesize(scheme);
    
    rmse(n) = norm(s1-s2)/sqrt(nScheme);
    
    % correct parameters
    p = model.getParams('excludeS0', true);
    testData(n,1).comp1_vFraction = p(1);
    testData(n,1).comp1_diff      = p(2);
    
    testData(n,1).comp2_vFraction = 1-p(1);
    testData(n,1).comp2_diff      = p(3);
    testData(n,1).comp2_theta     = p(4);
    testData(n,1).comp2_phi       = p(5);
    
    
end            

if max(rmse) < tol
    fprintf('rotateAxis function works fine in class %s (max rmse = %1.2e).\n', model.name, max(rmse));
    save('/Users/medmfarb/Documents/mycodes/diffusionMRI/caminoInMatlab/test/testData/synthetic/ball_stick_rotated.mat', 'testData', 'scheme');
else
     error('rotateAxis function is buggy in class %s! (max rmse = %1.2e)\n', model.name, max(rmse));
end

%% cylinder_zeppelin
load(fullfile(path_to_data, 'cylinder_zeppelin.mat'), 'testData', 'scheme');
N = length(testData);
nScheme = size(scheme, 1);

rmse = zeros(N, 1);
for n = 1:100
    s0 = 1;
    
    % read the input parameters
    comp1 = cylinder([s0, testData(n).comp1_diff, testData(n).comp1_r, ...
                     testData(n).comp1_theta, testData(n).comp1_phi]);
               
    comp2 = zeppelin([s0, testData(n).comp2_diffPar, ...
                      testData(n).comp2_diffPerp, testData(n).comp2_theta,...
                      testData(n).comp2_phi]);
    
    
    vFractions = [testData(n).comp1_vFraction, testData(n).comp2_vFraction];
    model = twoCompartmentModel(comp1, comp2, 'vFractions', vFractions);  
    s1 = model.synthesize(scheme);
    
    % rotate the axis and synthesize signals again
    model.rotateAxis();
    s2 = model.synthesize(scheme);
    
    rmse(n) = norm(s1-s2)/sqrt(nScheme);
    
    % correct parameters
    p = model.getParams('excludeS0', true);
    testData(n,1).comp1_vFraction = p(1);
    testData(n,1).comp1_diff      = p(2);
    testData(n,1).comp1_r         = p(3);
    testData(n,1).comp1_theta     = p(4);
    testData(n,1).comp1_phi       = p(5);
    
    testData(n,1).comp2_vFraction = 1-p(1);
    testData(n,1).comp2_diffPar   = p(6);
    testData(n,1).comp2_diffPerp  = p(7);
    testData(n,1).comp2_theta     = p(8);
    testData(n,1).comp2_phi       = p(9);
    
end            

if max(rmse) < tol
    fprintf('rotateAxis function works fine in class %s (max rmse = %1.2e).\n', model.name, max(rmse));
    save('/Users/medmfarb/Documents/mycodes/diffusionMRI/caminoInMatlab/test/testData/synthetic/cylinder_zeppelin_rotated.mat', 'testData', 'scheme');
else
     error('rotateAxis function is buggy in class %s! (max rmse = %1.2e)\n', model.name, max(rmse));
end

%% ellipticalCylinder_ball_stick
load(fullfile(path_to_data, 'ellipticalCylinder_ball_stick.mat'), 'testData', 'scheme');
N = length(testData);
nScheme = size(scheme, 1);

rmse = zeros(N, 1);
for n = 1:100
    s0 = 1;
    
    % read the input parameters
    comp1 = ellipticalCylinder([s0, testData(n).comp1_diff, ...
                 testData(n).comp1_r1, testData(n).comp1_r2 ...
                 testData(n).comp1_theta, testData(n).comp1_phi, ...
                 testData(n).comp1_alpha]);
               
    comp2 = ball([s0, testData(n).comp2_diff]);
    
    comp3 = stick([s0, testData(n).comp3_diff, testData(n).comp3_theta,...
                   testData(n).comp3_phi]);
    
    
    vFractions = [testData(n).comp1_vFraction, ...
                  testData(n).comp2_vFraction, testData(n).comp3_vFraction];
              
    model = threeCompartmentModel(comp1, comp2, comp3, ...
                                  'vFractions', vFractions);  
    s1 = model.synthesize(scheme);
    
    % rotate the axis and synthesize signals again
    model.rotateAxis();
    s2 = model.synthesize(scheme);
    
    rmse(n) = norm(s1-s2)/sqrt(nScheme);
    
    % correct parameters
    p = model.getParams('excludeS0', true);
    testData(n,1).comp1_vFraction = p(1);
    testData(n,1).comp1_diff      = p(3);
    testData(n,1).comp1_r1        = p(4);
    testData(n,1).comp1_r2        = p(5);
    testData(n,1).comp1_theta     = p(6);
    testData(n,1).comp1_phi       = p(7);
    testData(n,1).comp1_alpha     = p(8);
    
    testData(n,1).comp2_vFraction = p(2);
    testData(n,1).comp2_diff      = p(9);
    
    testData(n,1).comp3_vFraction = 1 - p(1) - p(2);
    testData(n,1).comp3_diff      = p(10);
    testData(n,1).comp3_theta     = p(11);
    testData(n,1).comp3_phi       = p(12);
    
end            

if max(rmse) < tol
    fprintf('rotateAxis function works fine in class %s (max rmse = %1.2e).\n', model.name, max(rmse));
    save('/Users/medmfarb/Documents/mycodes/diffusionMRI/caminoInMatlab/test/testData/synthetic/ellitpicalCylinder_ball_stick_rotated.mat', 'testData', 'scheme');
else
     error('rotateAxis function is buggy in class %s! (max rmse = %1.2e)\n', model.name, max(rmse));
end