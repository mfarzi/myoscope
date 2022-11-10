clc; clear all; close all;
% add path to the myoscope toolbox
addpath(genpath(pwd));

%% Model Construction
% construct a three-compartment model with an elliptical 
% cylinder, ball, and stick
model = multicompartment.str2model('[cylinderECS-ball-stick]');

% print model parameters' name
paramsName = model.getParamsName();
fprintf('Parameters for %s\n', model.name)
fprintf('%s, ', paramsName{:})
fprintf('\n')

% add constraints to make stick and cylinder parallel
% Note: orientation is represented using Euler's angles
model.addConstraint('comp3.theta=comp1.theta')
model.addConstraint('comp3.phi=comp1.phi')

% add constraints on parallel diffusivity in cylinder
model.addConstraint('comp1.diffPar>=0.8e-9')
model.addConstraint('comp1.diffPar<=1.2e-9')

% add constriant to set ball diffusivity equal to the stick diffusivity
model.addConstraint('comp3.diffPar=comp2.diff')

% add constraint on the ball diffusivity
model.addConstraint('comp2.diff>=2e-9')
model.addConstraint('comp2.diff<=2.5e-9')

% normalise the s0 signal to 1
model.addConstraint('s0=1')

%% create a stejeskal-tanner diffusion scheme file
ghat = [1               0             0           ;
        0.2619951531    0.4319920082  0.8629840349;
        0.5431406917   -0.4881264411 -0.6831769658;
       -0.814710829    -0.3858630429  0.4328463668;
        0.08798196555   0.1849620867 -0.9787993667;
       -0.003999318174 -0.9098448847 -0.4149292606;
        0.5822408064    0.8003310053  0.1430591672;
        0.08599286289   0.866928048  -0.4909592521;
       -0.678931092     0.1389858936 -0.7209268296;
        0.6932992232   -0.6983013821  0.1780768567];
nDirections = size(ghat, 1);
bval = [0 69 280 620 1100 1700 2500]; % s/mm^2
dt = [10 30 50];                % diffusion time in ms
delta = 2.5;                          % ms
te = 55;                              % echo time in ms

schemefile = scheme();
for DT = dt
    for B = bval
        schemefile.add(ghat, ones(nDirections, 1)*B*1e6, DT*1e-3, delta*1e-3, te*1e-3)
    end
end
%% generate a parameter vector 
% uncomment code below for a random vector
% params = model.randomInit();

% or set it manually
s0 = 1;
comp1_vol = 0.6;
comp1_diffPar = 0.9e-9;
comp1_r1 = 12e-6;
comp1_r2 = 8e-6;
comp1_theta = 0.2;
comp1_phi = 0;
comp1_alpha = 0;
comp2_vol = 0.2;
comp2_diff = 2e-9;
comp3_vol = 1-comp1_vol-comp2_vol;
comp3_diffPar=comp2_diff;
comp3_theta = comp1_theta;
comp3_phi = comp1_phi; 
params = [s0, comp1_vol, comp1_diffPar, comp1_r1, comp1_r2, comp1_theta, ...
          comp1_phi, comp1_alpha, comp2_vol, comp2_diff, comp3_vol, ...
          comp3_diffPar, comp3_theta, comp3_phi]';

%% generate synthetic diffusion MRI signal
sig = model.synthesize(params, schemefile);

% plot average orientation signals at dt=10, 30 and 50 ms
nBvals = length(bval);
x = bval';

y10 = zeros(nBvals, 1);
for i=1:nBvals
    idx = schemefile.dt==0.01 & schemefile.bval==bval(i)*1e6;
    y10(i) = mean(sig(idx));
end

y30 = zeros(nBvals, 1);
for i=1:nBvals
    idx = schemefile.dt==0.03 & schemefile.bval==bval(i)*1e6;
    y30(i) = mean(sig(idx));
end

y50 = zeros(nBvals, 1);
for i=1:nBvals
    idx = schemefile.dt==0.05 & schemefile.bval==bval(i)*1e6;
    y50(i) = mean(sig(idx));
end
figure('position', [100 100 700 400]); 
hold on;
h1 = plot(x, y10, 'k','LineWidth', 3);
h2 = plot(x, y30, 'r','LineWidth', 3);
h3 = plot(x, y50, 'b','LineWidth', 3);
xlabel('b-value [s/mm^2]')
ylabel('Signal Attenuation')
legend([h1, h2, h3], {'DT=10ms', 'DT=30ms', 'DT=50ms'})
set(gca, 'fontsize', 18)


%% estimate model parameters from the input data
[p, rmse] = model.fit(sig, schemefile, 'nreps', 5);
fprintf('\nThe RMSE for paramter estimation is %f\n', rmse)

% print estimated model parameters
fprintf('\nEstimated model parameers versus Ground-Truth (SNR=Infinity):\n')
for i = 1: model.getParamsNum
    fprintf('%s: %1.4e versus %1.4e\n', model.getParamsName{i}, p(i), params(i));
end

%% estimate model parameters at the presence of Rician noise
noise_std = 0.03;
M = size(sig, 1);
sig_noisy = abs((sig+randn(M,1)*noise_std)+randn(M,1)*noise_std*1i);
[p, rmse] = model.fit(sig_noisy, schemefile, 'nreps', 5);
fprintf('\nThe RMSE for paramter estimation is %f\n', rmse)

% print estimated model parameters
fprintf('\nEstimated model parameers versus Ground-Truth (SNR=%2d dB):\n', floor(-20*log10(noise_std)))
for i = 1: model.getParamsNum
    fprintf('%s: %1.4e versus %1.4e\n', model.getParamsName{i}, p(i), params(i));
end
