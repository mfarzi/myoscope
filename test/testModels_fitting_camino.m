% This script test if the optimisation procedure converges to the given
% parameters for each model.

clc; clear all; close all
addpath('../models');
addpath('../');
addpath('/Users/medmfarb/Documents/mycodes/utility/');
addpath('../../../camino/matlab_help_functions');

path_to_data = ['/Users/medmfarb/Documents/mycodes/diffusionMRI/',...
                'caminoInMatlab/test/testData/synthetic/'];

tol = 1e-10;            

%% ellipticalCylinder with camino

% THE CURRENT IMPLEMENTATION IN CAMINO REQUIRES A TWO PHASE OPTIMIZATION.
% FIRST, A DIFFUSION TENSOR IS FITTED AND THEN THE SCHEME FILE IS ROTATED.
% NEXT, PARAMETERS THETA, PHI, AND ALPHA ARE FIXED AND THE REST ARE
% ESTIMATED. 

% IN CAMINO IMPELMENTATION, DIFFUSIVITY OF BALL AND STICK ARE ALSO FIXED.
% SO DIRECT VALIDATION USING THE SAME DATASET USED FOR MATLAB
% IMPLEMENTATION IS NOT POSSIBLE.

modelName = 'ellipticalCylinder';
nParams = 6;

load(fullfile(path_to_data, strcat(modelName, '_rotated.mat')), 'testData', 'scheme');
camino_write_scheme(scheme(:,1:7), './tmp/scheme.scheme');

nScheme = size(scheme, 1);
nMultiRun = 5;
N = length(testData);
err = zeros(N, nParams);
rmse = zeros(N, 1);

P1 = zeros(N, nParams);
P2 = zeros(N, nParams);

hbar = parfor_progressbar(N, ...
            sprintf('Optimisation is running for %s ...', modelName));
tic;        
parfor n = 1:N
    % read the input parameters
    p1 = [testData(n).diff, testData(n).r1, testData(n).r2,...
          testData(n).theta, testData(n).phi, testData(n).alpha];
      
    P1(n,:) = p1;       
    
    % write data for camino
    filename = sprintf('./tmp/data%d.Bdouble', n);
    fid = fopen(filename, 'w', 'b');
    fwrite(fid, testData(n).s, 'double');
    fclose(fid);
    
    % model fitting from random initial condition
    str = ['modelfit -inputfile ', filename, ' -inputdatatype double ',...
           '-fitalgorithm MULTIRUNLM -samples ', num2str(nMultiRun), ' ', ... % '-fitalgorithm LM ', ...
           '-fitmodel ELLIPTICALCYLINDERONLY ', ...
           '-schemefile ./tmp/scheme.scheme ',...
           '-noisemodel Gaussian ', ...
           '-voxels 1 ',...
           '-outputfile ', sprintf('./tmp/params%d.Bdouble', n)];
    system(str);

    % read the input
    fid= fopen(sprintf('./tmp/params%d.Bdouble', n),'r','b');
    d=fread(fid,'double');
    fclose(fid);
    d = reshape(d, [length(d)/nMultiRun, nMultiRun]);
    [~,iMIN] = min(d(end,:));
    d = d(:,iMIN); 
    
    p2 = [d(4), d(7), d(8), d(5), d(6), d(9)];
    model = ellipticalCylinder([1, p2]);
    model.rotateAxis();
    p2 = model.getParams('excludeS0', true);
    thisCost = model.getCost(scheme, testData(n).s);
    
    P2(n,:) = p2;
    rmse(n) = sqrt(thisCost/nScheme);
    
    err(n, :) = [abs(p2(1)-p1(1))/p1(1) * 100, ...
                 abs(p2(2)-p1(2))/p1(2) * 100, ...
                 abs(p2(3)-p1(3))/p1(3) * 100, ...
                 min(mod(p2(4)-p1(4), pi), mod(p1(4)-p2(4), pi)), ...
                 min(mod(p2(5)-p1(5), 2*pi), mod(p1(5)-p2(5), 2*pi)),...
                 min(mod(p2(6)-p1(6), pi), mod(p1(6)-p2(6), pi))];
     
    hbar.iterate(1); % update progress by one iteration
end
time = toc;
close(hbar); % close the progress bar

%%%%%%%%%%%%%%%%%%%%%%    Convergence Criteria    %%%%%%%%%%%%%%%%%%%%%%%%%
% for diffusivity and radius, the normaised difference should be less     %
% than 1 percent.                                                         %
% for angles, the difference should be less than 10 degree (0.1745 radian)%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
isConverged = rmse<1e-4 & ...
              err(:,1)<1 & err(:,2)<1 & err(:,3)<1 & ...
              err(:,4) < 0.01745 & err(:,5) < 0.01745 & err(:,6) < 0.01745;

if all(isConverged)
    fprintf('optimisation module works fine in class %s. (time per case = %1.1f milisecond)\n', modelName, time/N*1000);
else
    fprintf('optimissation converged for %d/%d tests for %s model. (time per case = %1.1f milisecond)\n', sum(isConverged), N, modelName, time/N*1000);
end

%% ellipticalCylinderGPD_ball_stick with camino
% THE EXISTING CODE CANNOT WORK
% modelName = 'ellipticalCylinder_ball_stick';
% nParams = 12;
% 
% load(fullfile(path_to_data, strcat(modelName, '_rotated.mat')), 'testData', 'scheme');
% nScheme = size(scheme, 1);
% nMultiRun = 10;
% N = length(testData);
% err = ones(N, nParams);
% rmse = ones(N, 1);
% flag = zeros(N, 1);
% 
% P1 = zeros(N, nParams);
% P2 = zeros(N, nParams);
% 
% hbar = parfor_progressbar(N, ...
%             sprintf('Optimisation is running for %s ...', modelName));
% tic;        
% for n = 33%1:N
%     % read the input parameters
%     s0 = 1;
%     p1a = [testData(n).comp1_diff, testData(n).comp1_r1   , ...
%            testData(n).comp1_r2  , testData(n).comp1_theta, ...
%            testData(n).comp1_phi , testData(n).comp1_alpha];
%     p1b =  testData(n).comp2_diff;
%     p1c = [testData(n).comp3_diff, testData(n).comp3_theta, ...
%            testData(n).comp3_phi];
%     vFractions = [testData(n).comp1_vFraction, ...
%                   testData(n).comp2_vFraction, testData(n).comp3_vFraction]; 
%     
%     p1 = [vFractions(1:2), p1a, p1b, p1c]';
%     P1(n,:) = p1;  
%     
%     % write data for camino
%     filename = sprintf('./tmp/data%d.Bdouble', n);
%     fid = fopen(filename, 'w', 'b');
%     fwrite(fid, testData(n).s, 'double');
%     fclose(fid);
%      
%     % model fitting from random initial condition
%     str = ['modelfit -inputfile ', filename, ' -inputdatatype double ',...
%            '-fitalgorithm MULTIRUNLM -samples ', num2str(nMultiRun), ' ', ... % '-fitalgorithm LM ', ...
%            '-fitmodel ELLIPTICALCYLINDERTENSORSTICK ', ...
%            '-schemefile ./tmp/scheme.scheme ',...
%            '-noisemodel Gaussian ', ...
%            '-voxels 1 ',...
%            '-outputfile ', sprintf('./tmp/params%d.Bdouble', n)];
%     system(str);
% 
%     % read the input
%     fid= fopen(sprintf('./tmp/params%d.Bdouble', n),'r','b');
%     d=fread(fid,'double');
%     fclose(fid);
%     d = reshape(d, [length(d)/nMultiRun, nMultiRun]);
%     [~,iMIN] = min(d(end,:));
%     d = d(:,iMIN); 
%     
%     p2 = [d(3), d(4), d(6), d(9), d(10), d(7), d(8), d(11), d(12), d(18), d(19), d(20)];
%     p2a = p2(3:8);
%     p2b = p2(9);
%     p2c = p2(10:12);
%     
%     comp1 = ellipticalCylinder([s0, p2a]);
%     comp2 = ball([s0, p2b]);
%     comp3 = stick([s0, p2c]);
%     
%     % model fitting from random initial condition
%     model = threeCompartmentModel(comp1, comp2, comp3, 'vFractions', vFractions);
%     model.rotateAxis();
%     
%     p2 = model.getParams('excludeS0', true);
%     P2(n,:) = p2;
% 
%     
%     flag(n) = d(1);
%     rmse(n) = sqrt(model.getCost(scheme, testData(n).s)/nScheme);
%     
%     err(n, :) = abs(p1-p2)./p1*100;
%      
%     hbar.iterate(1); % update progress by one iteration
% end
% time = toc;
% close(hbar); % close the progress bar
% 
% %%%%%%%%%%%%%%%%%%%%%%    Convergence Criteria    %%%%%%%%%%%%%%%%%%%%%%%%%
% % for diffusivity and radius, the normaised difference should be less     %
% % than 1 percent.                                                         %
% % for angles, the difference should be less than 10 degree (0.1745 radian)%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% isConverged = rmse<1e-4 & sum(err<1, 2)==nParams;
% 
% if all(isConverged)
%     fprintf('optimisation module works fine in class %s. (time per case = %1.1f milisecond)\n', modelName, time/N*1000);
% else
%     fprintf('optimissation converged for %d/%d tests for %s model. (time per case = %1.1f milisecond)\n', sum(isConverged), N, modelName, time/N*1000);
% end