clc; clear all; close all;
addpath('../models');
addpath('../');
addpath('../linker');
addpath('../utility');
addpath('../../utility/');

load('/Users/medmfarb/Documents/mycodes/diffusionMRI/caminoInMatlab/test/testData/experimental/scheme.mat', 'scheme');
nBins = 300;
s = zeros(305, 1);
for theta = linspace(0, pi, nBins)
    model = cylinderGPD('s0', 1, 'diffPar', 2e-9, 'r', 1e-6, 'theta', theta, 'phi', pi/4);
    
    s = s + model.synth(scheme);
end

s = s/nBins;

t = tensor();
t.addConstraint('s0=1');
t.addConstraint(sprintf('phi = %f', pi/4));
t.fitMultiRun(scheme, s, 10);
