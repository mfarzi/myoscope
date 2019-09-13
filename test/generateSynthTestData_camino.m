% This script generates synthetic data for different compartment models
% with various initial parameters using camino. The results will be used
% as ground truth for validation of the developed matlab models.
clc; clear all; close all;



addpath('../../../camino/matlab_help_functions');
addpath('../models');

load('/Users/medmfarb/Documents/mycodes/diffusionMRI/caminoInMatlab/test/testData/synthetic/scheme.mat');

delete('./tmp/*');
% write the scheme file into the tmp folder
camino_write_scheme(scheme(:,1:7), './tmp/scheme.scheme');

%% TENSOR
rng(1234);                              % set the random see
N = 100;                                % numober of total tests
for n = 1:100

    diff      = (sort(rand(3,1), 'descend')*2.5 + 0.5)*1e-9;
    diffPar   = diff(1);
    diffPerp1 = diff(2);
    diffPerp2 = diff(3);
    theta     = rand(1)*2*pi;
    phi       = rand(1)*2*pi;
    alpha     = rand(1)*2*pi;
    
    str = sprintf(['datasynth -synthmodel compartment 1 ' ...
                   'TENSOR %d %d %d %d %d %d ',...
                   '-schemefile ./tmp/scheme.scheme ',...
                   '-voxels 1 -outputfile ./tmp/synth_data.Bfloat'], ...
                    diffPar, theta, phi, diffPerp1, diffPerp2, alpha); % 
    system(str);

    % read the synthetic data
    fid= fopen('./tmp/synth_data.Bfloat','r','b');
    s =fread(fid,'float');
    fclose(fid);

    testData(n,1).diffPar   = diffPar;
    testData(n,1).diffPerp1 = diffPerp1;
    testData(n,1).diffPerp2 = diffPerp2;
    testData(n,1).theta     = theta;
    testData(n,1).phi       = phi;
    testData(n,1).alpha     = alpha;
    testData(n,1).s         = s;

    delete('./tmp/synth_data.Bfloat');
end
save('/Users/medmfarb/Documents/mycodes/diffusionMRI/caminoInMatlab/test/testData/synthetic/tensor.mat', 'testData', 'scheme');
clear testData

%% BALL
rng(1234);                              % set the random see
N = 100;                                % numober of total tests
for n = 1:N
    diff = (rand(1)*2.5+0.5)*1e-9;
    
    str = sprintf(['datasynth -synthmodel compartment 1 ' ...
                   'BALL %d ',...
                   '-schemefile ./tmp/scheme.scheme ',...
                   '-voxels 1 -outputfile ./tmp/synth_data.Bfloat'], ...
                    diff); % 
    system(str);

    % read the synthetic data
    fid= fopen('./tmp/synth_data.Bfloat','r','b');
    s =fread(fid,'float');
    fclose(fid);

    testData(n,1).diff = diff;
    testData(n,1).s    = s;
    
    delete('./tmp/synth_data.Bfloat');
end

save('/Users/medmfarb/Documents/mycodes/diffusionMRI/caminoInMatlab/test/testData/synthetic/ball.mat', 'testData', 'scheme');
clear testData

%% Zeppelin
rng(1234);                              % set the random see
N = 100;                                % numober of total tests
for n = 1:100
    diff     = (sort(rand(2,1), 'descend')*2.5 + 0.5)*1e-9;
    diffPar  = diff(1);
    diffPerp = diff(2);
    theta    = rand(1)*2*pi;
    phi      = rand(1)*2*pi;
    
    str      = sprintf(['datasynth -synthmodel compartment 1 ' ...
                   'ZEPPELIN %d %d %d %d ',...
                   '-schemefile ./tmp/scheme.scheme ',...
                   '-voxels 1 -outputfile ./tmp/synth_data.Bfloat'], ...
                    diffPar, theta, phi, diffPerp); % 
    system(str);

    % read the synthetic data
    fid= fopen('./tmp/synth_data.Bfloat','r','b');
    s =fread(fid,'float');
    fclose(fid);

    testData(n,1).diffPar  = diffPar;
    testData(n,1).diffPerp = diffPerp;
    testData(n,1).theta    = theta;
    testData(n,1).phi      = phi;
    testData(n,1).s        = s;
    
    delete('./tmp/synth_data.Bfloat');
end

save('/Users/medmfarb/Documents/mycodes/diffusionMRI/caminoInMatlab/test/testData/synthetic/zeppelin.mat', 'testData', 'scheme');
clear testData

%% STICK
rng(1234);                              % set the random see
N = 100;                                % numober of total tests
for n = 1:N
    diff = (rand(1)*2.5+0.5)*1e-9;
    theta = rand(1)*2*pi;
    phi = rand(1)*2*pi;
    
    str = sprintf(['datasynth -synthmodel compartment 1 ' ...
                   'STICK %d %d %d ',...
                   '-schemefile ./tmp/scheme.scheme ',...
                   '-voxels 1 -outputfile ./tmp/synth_data.Bfloat'], ...
                    diff, theta, phi); % 
    system(str);

    % read the synthetic data
    fid= fopen('./tmp/synth_data.Bfloat','r','b');
    s =fread(fid,'float');
    fclose(fid);

    testData(n,1).diff  = diff;
    testData(n,1).theta = theta;
    testData(n,1).phi   = phi;
    testData(n,1).s     = s;
    
    delete('./tmp/synth_data.Bfloat');
end

save('/Users/medmfarb/Documents/mycodes/diffusionMRI/caminoInMatlab/test/testData/synthetic/stick.mat', 'testData', 'scheme');
clear testData;

%% CYLINDER
rng(1234);                              % set the random see
N = 100;                                % numober of total tests
for n = 1:100
    
    diff  = (rand(1)* 2.5 + 0.5)*1e-9;  % random number in range [0.5e-9, 3e-9] m^2/s
    r     = (rand(1)*19.5 + 0.5)*1e-6;  % random number in range [0.5e-6, 20e-6] m
    theta = rand(1)*2*pi;
    phi   = rand(1)*2*pi;
    str   = sprintf(['datasynth -synthmodel compartment 1 ' ...
                     'CYLINDERGPD %d %d %d %d ',...
                     '-schemefile ./tmp/scheme.scheme ',...
                     '-voxels 1 -outputfile ./tmp/synth_data.Bfloat'], ...
                     diff, theta, phi, r); % 
    system(str);

    % read the synthetic data
    fid= fopen('./tmp/synth_data.Bfloat','r','b');
    s =fread(fid,'float');
    fclose(fid);

    testData(n,1).diff  = diff;
    testData(n,1).r     = r;
    testData(n,1).theta = theta;
    testData(n,1).phi   = phi;
    testData(n,1).s     = s;
    
    delete('./tmp/synth_data.Bfloat');
end

save('/Users/medmfarb/Documents/mycodes/diffusionMRI/caminoInMatlab/test/testData/synthetic/cylinder.mat', 'testData', 'scheme');
clear testData;

%% CYLINDERGDR
rng(1234);                              % set the random see
N = 100;                                % numober of total tests
for n = 1:100
    
    diff  = (rand(1)* 2.5 + 0.5)*1e-9;  % random number in range [0.5e-9, 3e-9] m^2/s
    r     = (rand(1)*19.5 + 0.5)*1e-6;       % random number in range [0.5, 20] um
    kappa = rand(1)*20;             % random kappa values
    nu    = r/kappa;
    theta = rand(1)*pi/2;
    phi   = rand(1)*2*pi;
    str   = sprintf(['datasynth -synthmodel compartment 1 ' ...
                     'GAMMADISTRIBRADIICYLINDERS %d %d %d %d %d ',...
                     '-schemefile ./tmp/scheme.scheme ',...
                     '-voxels 1 -outputfile ./tmp/synth_data.Bfloat'], ...
                     kappa, nu, diff, theta, phi); % 
    system(str);

    % read the synthetic data
    fid= fopen('./tmp/synth_data.Bfloat','r','b');
    s =fread(fid,'float');
    fclose(fid);

    testData(n,1).diff  = diff;
    testData(n,1).kappa = kappa;
    testData(n,1).nu    = nu;
    testData(n,1).theta = theta;
    testData(n,1).phi   = phi;
    testData(n,1).s     = s;
    
    delete('./tmp/synth_data.Bfloat');
end

save('/Users/medmfarb/Documents/mycodes/diffusionMRI/caminoInMatlab/test/testData/synthetic/cylinderGDR.mat', 'testData', 'scheme');
clear testData;

%% ELLIPTICALCYLINDER
rng(1234);                              % set the random see
N = 100;                                % numober of total tests
for n = 1:N
    diff  = (rand(1)* 2.5 + 0.5)*1e-9;  % random number in range [0.5e-9, 3e-9] m^2/s
    r     = (sort(rand(2,1), 'descend')*19.5 + 0.5)*1e-6;
    r1    = r(1);                       % random number in range [0.5e-6, 20e-6] m
    r2    = r(2);                       % random number in range [0.5e-6, 20e-6] m
    theta = rand(1)*2*pi;
    phi   = rand(1)*2*pi;
    alpha = rand(1)*2*pi;
    
    str   = sprintf(['datasynth -synthmodel compartment 1 '           ...
                     'ELLIPTICALCYLINDER %d %d %d %d %d %d '         ,...
                     '-schemefile ./tmp/scheme.scheme '              ,...
                     '-voxels 1 -outputfile ./tmp/synth_data.Bfloat'],...
                     diff, theta, phi, r1, r2, alpha); % 
    system(str);

    % read the synthetic data
    fid= fopen('./tmp/synth_data.Bfloat','r','b');
    s =fread(fid,'float');
    fclose(fid);

    testData(n,1).diff  = diff;
    testData(n,1).r1    = r1;
    testData(n,1).r2    = r2;
    testData(n,1).theta = theta;
    testData(n,1).phi   = phi;
    testData(n,1).alpha = alpha;
    testData(n,1).s     = s;
    
    delete('./tmp/synth_data.Bfloat');
end

save('/Users/medmfarb/Documents/mycodes/diffusionMRI/caminoInMatlab/test/testData/synthetic/ellipticalCylinder.mat', 'testData', 'scheme');
clear testData;

%% BALL_STICK
rng(1234);                              % set the random see
N = 100;                                % numober of total tests
for n = 1:N
    v = rand(2,1); v = v/sum(v);
    
    v_ball      = v(1);
    v_stick     = v(2);
    
    diff_ball   = (rand(1)*2.5+0.5)*1e-9;
    
    diff_stick  = (rand(1)*2.5+0.5)*1e-9;
    theta_stick = rand(1)*2*pi;
    phi_stick   = rand(1)*2*pi;
    
    str = sprintf(['datasynth -synthmodel compartment 2 ' ...
                   'BALL %d %d ', ...
                   'STICK %d %d %d ',...
                   '-schemefile ./tmp/scheme.scheme ',...
                   '-voxels 1 -outputfile ./tmp/synth_data.Bfloat'], ...
                    v_ball, diff_ball, ...
                    diff_stick, theta_stick, phi_stick); % 
    system(str);

    % read the synthetic data
    fid= fopen('./tmp/synth_data.Bfloat','r','b');
    s =fread(fid,'float');
    fclose(fid);

    testData(n,1).comp1_diff      = diff_ball;
    testData(n,1).comp1_vFraction = v_ball;
    
    testData(n,1).comp2_diff      = diff_stick;
    testData(n,1).comp2_theta     = theta_stick;
    testData(n,1).comp2_phi       = phi_stick;
    testData(n,1).comp2_vFraction = v_stick;
    
    testData(n,1).s = s;
    
    delete('./tmp/synth_data.Bfloat');
end

save('/Users/medmfarb/Documents/mycodes/diffusionMRI/caminoInMatlab/test/testData/synthetic/ball_stick.mat', 'testData', 'scheme');
clear testData;

%% CYLINDER_ZEPPELIN
rng(1234);                              % set the random see
N = 100;                                % numober of total tests
for n = 1:N
    
    v = rand(2,1); v = v/sum(v);
    
    v_cylinder      = v(1);
    v_zeppelin     = v(2);
    
    diff_cylinder  = (rand(1)*2.5+0.5)*1e-9;
    r_cylinder     = (rand(1)*19.5 + 0.5)*1e-6;
    theta_cylinder = rand(1)*2*pi;
    phi_cylinder   = rand(1)*2*pi;
    
    diff = (sort(rand(2,1), 'descend')*2.5 + 0.5)*1e-9;
    diffPar_zeppelin  = diff(1);
    diffPerp_zeppelin = diff(2);
    theta_zeppelin    = rand(1)*2*pi;
    phi_zeppelin      = rand(1)*2*pi;
    
    str = sprintf(['datasynth -synthmodel compartment 2 ' ...
                   'CYLINDERGPD %d %d %d %d %d ',...
                   'ZEPPELIN %d %d %d %d ', ...
                   '-schemefile ./tmp/scheme.scheme ',...
                   '-voxels 1 -outputfile ./tmp/synth_data.Bfloat'], ...
                    v_cylinder, diff_cylinder, theta_cylinder, phi_cylinder, r_cylinder, ....
                    diffPar_zeppelin, theta_zeppelin, phi_zeppelin, diffPerp_zeppelin);
    system(str);

    % read the synthetic data
    fid= fopen('./tmp/synth_data.Bfloat','r','b');
    s =fread(fid,'float');
    fclose(fid);

    testData(n,1).comp1_diff      = diff_cylinder;
    testData(n,1).comp1_r         = r_cylinder;
    testData(n,1).comp1_theta     = theta_cylinder;
    testData(n,1).comp1_phi       = phi_cylinder;
    testData(n,1).comp1_vFraction = v_cylinder;
    
    testData(n,1).comp2_diffPar   = diffPar_zeppelin;
    testData(n,1).comp2_diffPerp  = diffPerp_zeppelin;
    testData(n,1).comp2_theta     = theta_zeppelin;
    testData(n,1).comp2_phi       = phi_zeppelin;
    testData(n,1).comp2_vFraction = v_zeppelin;
    
    testData(n,1).s = s;
    
    delete('./tmp/synth_data.Bfloat');
end

save('/Users/medmfarb/Documents/mycodes/diffusionMRI/caminoInMatlab/test/testData/synthetic/cylinder_zeppelin.mat', 'testData', 'scheme');
clear testData

%% ELLIPTICALCYLINDER_BALL_STICK
rng(1234);                              % set the random see
N = 100;                                % numober of total tests
for n = 1:N
    
    v = rand(3,1); v = v/sum(v);
    v_cylinder = v(1);
    v_ball     = v(2);
    v_stick    = v(3);
    
    diff_cylinder  = (rand(1)*2.5+0.5)*1e-9;
    r = (sort(rand(2,1), 'descend')*19.5 + 0.5)*1e-6;
    r1_cylinder    = r(1);
    r2_cylinder    = r(2);
    theta_cylinder = rand(1)*2*pi;
    phi_cylinder   = rand(1)*2*pi;
    alpha_cylinder = rand(1)*2*pi;
    
    diff_ball = (rand(1)*2.5+0.5)*1e-9;
    
    diff_stick = (rand(1)*2.5+0.5)*1e-9;
    theta_stick = rand(1)*2*pi;
    phi_stick = rand(1)*2*pi;
    
    str = sprintf(['datasynth -synthmodel compartment 3 ' ...
                   'ELLIPTICALCYLINDER %d %d %d %d %d %d %d ',...
                   'BALL %d %d ', ...
                   'STICK %d %d %d ',...
                   '-schemefile ./tmp/scheme.scheme ',...
                   '-voxels 1 -outputfile ./tmp/synth_data.Bfloat'], ...
                    v_cylinder, diff_cylinder, theta_cylinder, phi_cylinder, r1_cylinder, r2_cylinder, alpha_cylinder, ....
                    v_ball, diff_ball, ...
                    diff_stick, theta_stick, phi_stick); % 
    system(str);

    % read the synthetic data
    fid= fopen('./tmp/synth_data.Bfloat','r','b');
    s =fread(fid,'float');
    fclose(fid);

    testData(n,1).comp1_diff      = diff_cylinder;
    testData(n,1).comp1_r1        = r1_cylinder;
    testData(n,1).comp1_r2        = r2_cylinder;
    testData(n,1).comp1_theta     = theta_cylinder;
    testData(n,1).comp1_phi       = phi_cylinder;
    testData(n,1).comp1_alpha     = alpha_cylinder;
    testData(n,1).comp1_vFraction = v_cylinder;
    
    testData(n,1).comp2_diff      = diff_ball;
    testData(n,1).comp2_vFraction = v_ball;
    
    testData(n,1).comp3_diff      = diff_stick;
    testData(n,1).comp3_theta     = theta_stick;
    testData(n,1).comp3_phi       = phi_stick;
    testData(n,1).comp3_vFraction = v_stick;
    testData(n,1).s = s;
    
    delete('./tmp/synth_data.Bfloat');
end
save('/Users/medmfarb/Documents/mycodes/diffusionMRI/caminoInMatlab/test/testData/synthetic/ellipticalCylinder_ball_stick.mat', 'testData', 'scheme');
clear testData;

