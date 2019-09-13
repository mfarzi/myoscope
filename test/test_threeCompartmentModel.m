% this script test the ability of camino java implementation to estimate
% model parameters for a three compartment model:
% ellipticalCylinder-ball-stick
clc; clear all; close all;
addpath('../../../camino/matlab_help_functions');
addpath('../models');
addpath('../');

load('/Users/medmfarb/Documents/mycodes/diffusionMRI/caminoInMatlab/test/testData/synthetic/scheme.mat');
%rng(77);
N = 100;
P1 = zeros(6, N);
P2 = zeros(6, N);
convergence = zeros(2,N);
for i = 1
    %% write scheme file
    delete('./tmp/*');
    % write the scheme file into the tmp folder
    camino_write_scheme(scheme(:,1:7), './tmp/scheme.scheme');

    %% synthesize signal
    diff_cylinder = (rand(1)+0.5)*1e-9;
    r1_cylinder = (rand(1)*19.5+0.5)*1e-6;
    r2_cylinder = max(rand(1)*r1_cylinder, 0.5e-6);
    theta_cylinder = pi/2*rand(1);
    phi_cylinder = rand(1)*2*pi;
    alpha_cylinder = pi*rand(1);
    v_cylinder = 0.1;%rand(1);
    
    diff_stick = 1.6e-9;
    theta_stick = theta_cylinder;
    phi_stick = phi_cylinder;
    v_stick = min([rand(1)*(1-v_cylinder), v_cylinder, 0.2]);

    diff_ball = 1.6e-9;
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
    
    % read synth data
    fid= fopen('./tmp/data.Bfloat','r','b');
    synth_data=fread(fid,'float');
    fclose(fid); 
    
    % write out the input data
    which_data = scheme.DELTA == 0.01;
    filename = fullfile('./tmp/data_DT.Bfloat');
    fid = fopen(filename, 'w', 'b');
    fwrite(fid, synth_data(which_data), 'float');
    fclose(fid);
      
    % write the scheme file to be used by camino
    path2scheme = fullfile('./tmp/scheme_fit_DT.scheme');
    camino_write_scheme(scheme(which_data,1:7), path2scheme);
    
    % fit the diffusion tensor to the data using modelfit
    str = ['modelfit -inputfile ./tmp/data_DT.Bfloat -schemefile ./tmp/scheme_fit_DT.scheme',...
           ' -model ldt -outputfile ./tmp/DT.Bdouble'];
    system(str);
    
    % rotate the scheme file around the diffusion tensor
    fid= fopen('./tmp/DT.Bdouble','r','b');
    d=fread(fid,'double');
    fclose(fid);
    DT = [[d(3), d(4), d(5)];...
          [d(4), d(6), d(7)];...
          [d(5), d(7), d(8)]];
    
    [U, E] = eig(DT);
    E = diag(E);
    [~, idx] = sort(E, 'descend');
    U = U(:, idx);
    E = E(idx);  
    

    thisScheme = scheme(:,1:7);
    rotated_dir = [scheme.x, scheme.y, scheme.z]*U;
    thisScheme.x = rotated_dir(:,1);
    thisScheme.y = rotated_dir(:,2);
    thisScheme.z = rotated_dir(:,3);
    
    camino_write_scheme(thisScheme(:,1:7), './tmp/thisScheme.scheme');

    p1 = [v_cylinder; v_ball; v_stick; diff_cylinder; r1_cylinder; r2_cylinder];
    P1(:,i) = p1;
    %
    %% estimate parameters
    nReps = 10;
    str = ['modelfit -inputfile ./tmp/data.Bfloat -inputdatatype float ',...
           '-fitalgorithm MULTIRUNLM -samples ', num2str(nReps), ' ', ... %'-fitalgorithm LM ', ...
           '-fitmodel ELLIPTICALCYLINDERTENSORSTICK ', ...
           '-schemefile ./tmp/thisScheme.scheme ',...
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