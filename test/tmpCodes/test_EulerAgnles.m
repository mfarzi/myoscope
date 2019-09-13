% This script aims to test if the functins getUnitFrame and getEulerAngles
% works properly or not.
clc; clear all; close all;

theta = rand(1)*pi;
phi = rand(1)*2*pi - pi;
alpha = rand(1)*2*pi - pi;

U = getUnitFrame(theta, phi, alpha);

[t, p, a] = getEulerAngles(U);