% this script test the ability of camino java implementation to estimate
% model parameters for a three compartment model:
% ellipticalCylinder-ball-stick
clc; clear all; close all;
addpath('../../../camino/matlab_help_functions');
addpath('../models');
addpath('../');

load('/Users/medmfarb/Documents/mycodes/diffusionMRI/caminoInMatlab/test/testData/synthetic/scheme.mat');
rng(77);
N = 100;
P1 = zeros(6, N);
P2 = zeros(6, N);
convergence = zeros(2,N);
for i = 1:N
    %% write scheme file
    delete('./tmp/*');
    % write the scheme file into the tmp folder
    camino_write_scheme(scheme(:,1:7), './tmp/scheme.scheme');

    %% synthesize signal
    diff_cylinder = (rand(1)+0.5)*1e-9;
    r1_cylinder = (rand(1)*20+0.5)*1e-6;
    r2_cylinder = max(rand(1)*r1_cylinder, 0.5e-6);
    theta_cylinder = pi/2+rand(1)*0.2;
    phi_cylinder = 0+rand(1)*0.2;
    alpha_cylinder = pi/2+rand(1)*0.2;
    v_cylinder = rand(1);
    
    diff_stick = 2.3e-9;
    theta_stick = theta_cylinder;
    phi_stick = phi_cylinder;
    v_stick = min([rand(1)*(1-v_cylinder), v_cylinder, 0.2]);

    diff_ball = 2.3e-9;
    v_ball = 1-v_stick-v_cylinder;

    

    str = sprintf(['datasynth -synthmodel compartment 3 ' ...
                   'ELLIPTICALCYLINDER %d %d %d %d %d %d %d ',...
                   'BALL %d %d ', ...
                   'STICK %d %d %d ',...
                   '-schemefile ./tmp/scheme.scheme ',...
                   '-voxels 1 -outputfile ./tmp/data.Bfloat'], ...
                    v_cylinder, diff_cylinder, theta_cylinder, phi_cylinder, r1_cylinder, r2_cylinder, alpha_cylinder, ....
                    v_ball, diff_ball, ...
                    diff_stick, theta_stick, phi_stick); % 
    system(str);

%     comp1 = ellipticalCylinder([diff_cylinder, r1_cylinder, r2_cylinder, theta_cylinder, phi_cylinder, alpha_cylinder]);
%     comp2 = ball(diff_ball);
%     comp3 = stick([diff_stick, theta_stick, phi_stick]);
%     vFractions = [v_cylinder, v_ball, v_stick]';
%     model = threeCompartmentModel(comp1, comp2, comp3, 'vFractions', vFractions);
%     s = model.synthesize(scheme);
%     fid = fopen('./tmp/data.Bfloat', 'w', 'b');
%     fwrite(fid, s, 'float');
%     fclose(fid);
    p1 = [v_cylinder; v_ball; v_stick; diff_cylinder; r1_cylinder; r2_cylinder];
    P1(:,i) = p1;
    %
    %% estimate parameters
    nReps = 1;
    str = ['modelfit -inputfile ./tmp/data.Bfloat -inputdatatype float ',...
           '-fitalgorithm LM ', ... %'-fitalgorithm MULTIRUNLM -samples ', num2str(nReps), ' ', ... 
           '-fitmodel ELLIPTICALCYLINDERTENSORSTICK ', ...
           '-schemefile ./tmp/scheme.scheme ',...
           '-noisemodel Gaussian ', ...
           '-voxels 1 ',...
           '-outputfile ./tmp/data_params.Bdouble'];
    system(str);

    fid= fopen('./tmp/data_params.Bdouble','r','b');
    d=fread(fid,'double');
    fclose(fid);

    D = reshape(d, [length(d)/nReps, nReps]);
    if sum(D(1,:)==0) == 0
        [~,iMIN] = min(D(end,:));
    else
        D(:,D(1,:)~=0)=[];
        [~,iMIN] = min(D(end,:));
    end
    params_camino = D(:,iMIN);

    p2 = params_camino([3:5,6,9,10]);
    P2(:,i) = p2;
    convergence(:,i) = [params_camino(1); params_camino(end)];
end