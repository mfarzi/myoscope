clc; clear all; close all;

addpath('../models');

f = @(x) x^2+2e9*x+1;
gf = @(x) 2*x+2e9;

x0 = 1;
optParam.alphaMax = 1;
optParam.c1 = 1e-4;
optParam.c2 = 0.9;
optParam.StopTol = 1e-8;
optParam.display = 1;

[x_min, f_min, iter] = ConjugateGradient(f, gf, x0, optParam);

%GSSTol = 1e-3;
%[x_min, f_min, iter] = SD_GSS(f, gf, x0, optParam.StopTol, GSSTol);