% This script test if the optimisation procedure converges to the given
% parameters for each model.

clc; clear all; close all
addpath('../models');
addpath('../');
addpath('../linker');
addpath('../utility');
addpath('../../utility/');

path_to_data = ['./testData/synthetic/'];

algorithm = 'conjugate-gradient';%'levenberg-marquardt'; %'interior-point' %;            
%% tensor  
modelName = 'tensor';
nParams   = 6;
nMultiRun = 10;

load(fullfile(path_to_data, strcat(modelName, '_rotated.mat')), 'testData', 'scheme');

nScheme = size(scheme, 1);
N    = length(testData);

flag = zeros(N, 1);
rmse = zeros(N, 1);
err  = zeros(N, nParams);

P1   = zeros(N, nParams);
P2   = zeros(N, nParams);

hbar = parfor_progressbar(N, ...
            sprintf('Optimisation is running for %s ...', modelName));
tic;        
parfor n = 1:100
    % read the input parameters
    p1 = [testData(n).diffPar, testData(n).diffPerp1, testData(n).diffPerp2, ...
          testData(n).theta, testData(n).phi, testData(n).alpha]';    
    P1(n,:) = p1;       
    
    % model fitting from random initial condition
    model = tensor('fitter', optimizer('algorithm', algorithm));
    model.addConstraint('s0=1');
    %model.addConstraint('diffPar>=diffPerp1');
    %model.addConstraint('diffPerp1>=diffPerp2');
    
    %model.randomInit();
    %p = model.fit(scheme, testData(n).s);
    
    p = model.fitMultiRun(scheme, testData(n).s, nMultiRun);
    
    flag(n) = p(1);
    thisCost = p(2);
    rmse(n) = sqrt(thisCost/nScheme);
    
    p2 = p(4:end); 
    P2(n,:) = p2;
  
    err(n, :) = abs(p1-p2)./p1*100;
     
    hbar.iterate(1); % update progress by one iteration
end
time = toc;
close(hbar); % close the progress bar

%%%%%%%%%%%%%%%%%%%%%%    Convergence Criteria    %%%%%%%%%%%%%%%%%%%%%%%%%
% a) The RMSE should be less than tol.                                    %
% b) Variation in model parameters should be less one percent.            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
isConverged = rmse<1e-4 & all(err<1, 2);

if all(isConverged)
    fprintf('optimisation module works fine in class %s. (time per case = %1.1f milisecond)\n', modelName, time/N*1000);
else
    fprintf('optimissation converged for %d/%d tests for %s model. (time per case = %1.1f milisecond)\n', sum(isConverged), N, modelName, time/N*1000);
end

%% ball   
modelName = 'ball';
nParams = 1;
nMultiRun = 5;

load(fullfile(path_to_data, strcat(modelName, '_rotated.mat')), 'testData', 'scheme');

nScheme = size(scheme, 1);
N    = length(testData);

flag = zeros(N, 1);
rmse = zeros(N, 1);
err  = zeros(N, nParams);

P1   = zeros(N, nParams);
P2   = zeros(N, nParams);

hbar = parfor_progressbar(N, ...
            sprintf('Optimisation is running for %s ...', modelName));
tic;        
parfor n = 1:N
    % read the input parameters
    p1 = testData(n).diff;
    P1(n,:) = p1;       
    
    % model fitting from random initial condition
    model = ball('fitter', optimizer('algorithm', algorithm));
    model.addConstraint('s0=1');
    
    %model.randomInit();
    %p = model.fit(scheme, testData(n).s);
    
    p = model.fitMultiRun(scheme, testData(n).s, nMultiRun);
    
    flag(n) = p(1);
    thisCost = p(2);
    rmse(n) = sqrt(thisCost/nScheme);
    
    p2 = p(4:end); 
    P2(n,:) = p2;
  
    err(n, :) = abs(p1-p2)./p1*100;
     
    hbar.iterate(1); % update progress by one iteration
end
time = toc;
close(hbar); % close the progress bar

%%%%%%%%%%%%%%%%%%%%%%    Convergence Criteria    %%%%%%%%%%%%%%%%%%%%%%%%%
% a) The RMSE should be less than tol.                                    %
% b) Variation in model parameters should be less one percent.            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
isConverged = rmse<1e-4 & all(err<1, 2);

if all(isConverged)
    fprintf('optimisation module works fine in class %s. (time per case = %1.1f milisecond)\n', modelName, time/N*1000);
else
    fprintf('optimissation converged for %d/%d tests for %s model. (time per case = %1.1f milisecond)\n', sum(isConverged), N, modelName, time/N*1000);
end

%% zeppelin
modelName = 'zeppelin';
nParams   = 4;
nMultiRun = 10;

load(fullfile(path_to_data, strcat(modelName, '_rotated.mat')), 'testData', 'scheme');

nScheme = size(scheme, 1);
N    = length(testData);

flag = zeros(N, 1);
rmse = zeros(N, 1);
err  = zeros(N, nParams);

P1   = zeros(N, nParams);
P2   = zeros(N, nParams);

hbar = parfor_progressbar(N, ...
            sprintf('Optimisation is running for %s ...', modelName));
tic;        
parfor n = 1:100
    % read the input parameters
    p1 = [testData(n).diffPar, testData(n).diffPerp, ...
          testData(n).theta, testData(n).phi]';    
    P1(n,:) = p1;       
    
    % model fitting from random initial condition
    model = zeppelin('s0', 1, 'fitter', optimizer('algorithm', algorithm));
    model.addConstraint('s0=1');
    
    %model.randomInit();
    %p = model.fit(scheme, testData(n).s);
    
    p = model.fitMultiRun(scheme, testData(n).s, nMultiRun);
    
    flag(n) = p(1);
    thisCost = p(2);
    rmse(n) = sqrt(thisCost/nScheme);
    
    p2 = p(4:end); 
    P2(n,:) = p2;
  
    err(n, :) = abs(p1-p2)./p1*100;
     
    hbar.iterate(1); % update progress by one iteration
end
time = toc;
close(hbar); % close the progress bar

%%%%%%%%%%%%%%%%%%%%%%    Convergence Criteria    %%%%%%%%%%%%%%%%%%%%%%%%%
% a) The RMSE should be less than tol.                                    %
% b) Variation in model parameters should be less one percent.            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
isConverged = rmse<1e-4 & all(err<1, 2);

if all(isConverged)
    fprintf('optimisation module works fine in class %s. (time per case = %1.1f milisecond)\n', modelName, time/N*1000);
else
    fprintf('optimissation converged for %d/%d tests for %s model. (time per case = %1.1f milisecond)\n', sum(isConverged), N, modelName, time/N*1000);
end

%% stick
modelName = 'stick';
nParams = 3;
nMultiRun = 10;

load(fullfile(path_to_data, strcat(modelName, '_rotated.mat')), 'testData', 'scheme');

nScheme = size(scheme, 1);
N    = length(testData);

flag = zeros(N, 1);
rmse = zeros(N, 1);
err  = zeros(N, nParams);

P1   = zeros(N, nParams);
P2   = zeros(N, nParams);

hbar = parfor_progressbar(N, ...
            sprintf('Optimisation is running for %s ...', modelName));
tic;        
parfor n = 1:N
    % read the input parameters
    p1 = [testData(n).diff, testData(n).theta, testData(n).phi]';
      
    P1(n,:) = p1;       
    
    % model fitting from random initial condition
    model = stick('fitter',  optimizer('algorithm', algorithm));
    model.addConstraint('s0=1');
    
    %model.randomInit();
    %p = model.fit(scheme, testData(n).s);
    
    p = model.fitMultiRun(scheme, testData(n).s, nMultiRun);
    
    flag(n) = p(1);
    thisCost = p(2);
    rmse(n) = sqrt(thisCost/nScheme);
    
    p2 = p(4:end); 
    P2(n,:) = p2;
  
    err(n, :) = abs(p1-p2)./p1*100;
     
    hbar.iterate(1); % update progress by one iteration
end
time = toc;
close(hbar); % close the progress bar

%%%%%%%%%%%%%%%%%%%%%%    Convergence Criteria    %%%%%%%%%%%%%%%%%%%%%%%%%
% a) The RMSE should be less than tol.                                    %
% b) Variation in model parameters should be less one percent.            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
isConverged = rmse<1e-4 & all(err<1, 2);

if all(isConverged)
    fprintf('optimisation module works fine in class %s. (time per case = %1.1f milisecond)\n', modelName, time/N*1000);
else
    fprintf('optimissation converged for %d/%d tests for %s model. (time per case = %1.1f milisecond)\n', sum(isConverged), N, modelName, time/N*1000);
end

%% cylinder
modelName = 'cylinder';
nParams = 4;
nMultiRun = 10;

load(fullfile(path_to_data, strcat(modelName, '_rotated.mat')), 'testData', 'scheme');

nScheme = size(scheme, 1);
N    = length(testData);

flag = zeros(N, 1);
rmse = zeros(N, 1);
err  = zeros(N, nParams);

P1   = zeros(N, nParams);
P2   = zeros(N, nParams);

hbar = parfor_progressbar(N, ...
            sprintf('Optimisation is running for %s ...', modelName));
tic;        
parfor n = 1:N
    % read the input parameters
    p1 = [testData(n).diff, testData(n).r, testData(n).theta, ...
          testData(n).phi]';
      
    P1(n,:) = p1;       
    
    % model fitting from random initial condition
    model = cylinder('fitter',  optimizer('algorithm', algorithm));
    model.addConstraint('s0=1');             
    
%     model.randomInit();
%     p = model.fit(scheme, testData(n).s);
    
    p = model.fitMultiRun(scheme, testData(n).s, nMultiRun);
    
    flag(n) = p(1);
    thisCost = p(2);
    rmse(n) = sqrt(thisCost/nScheme);
    
    p2 = p(4:end); 
    P2(n,:) = p2;
  
    err(n, :) = abs(p1-p2)./p1*100;
     
    hbar.iterate(1); % update progress by one iteration
end
time = toc;
close(hbar); % close the progress bar

%%%%%%%%%%%%%%%%%%%%%%    Convergence Criteria    %%%%%%%%%%%%%%%%%%%%%%%%%
% a) The RMSE should be less than tol.                                    %
% b) Variation in model parameters should be less one percent.            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
isConverged = rmse<1e-4 & all(err<1, 2);

if all(isConverged)
    fprintf('optimisation module works fine in class %s. (time per case = %1.1f milisecond)\n', modelName, time/N*1000);
else
    fprintf('optimissation converged for %d/%d tests for %s model. (time per case = %1.1f milisecond)\n', sum(isConverged), N, modelName, time/N*1000);
end

%% cylinder with Gamma Distributed Radii
modelName = 'cylinderGDR';
nParams = 5;
nMultiRun = 10;

load(fullfile(path_to_data, strcat(modelName, '.mat')), 'testData', 'scheme');

nScheme = size(scheme, 1);
N    = length(testData);

flag = zeros(N, 1);
rmse = zeros(N, 1);
err  = zeros(N, nParams);

P1   = zeros(N, nParams);
P2   = zeros(N, nParams);

hbar = parfor_progressbar(N, ...
            sprintf('Optimisation is running for %s ...', modelName));
tic;        
for n = 84%1:N
    % read the input parameters
    p1 = [1, testData(n).diff, testData(n).kappa, testData(n).nu, ...
          testData(n).theta, testData(n).phi]';
      
    P1(n,:) = p1(2:end);       
    
    % model fitting from random initial condition
    model = cylinderGDR(p1, 'fitter',  optimizer('algorithm', algorithm));
    model.addConstraint('s0=1');   
    thisSig = model.synth(scheme);
    
    model.addConstraint(sprintf('nu = %1.15f', p1(4)));   
    model.randomInit();
    
%     p = model.fit(scheme, thisSig);
    
    p = model.fitMultiRun(scheme, testData(n).s, nMultiRun);
    
    flag(n) = p(1);
    thisCost = p(2);
    rmse(n) = sqrt(thisCost/nScheme);
    
    p2 = p(4:end); 
    P2(n,:) = p2;
  
    err(n, :) = abs(p1(2:end)-p2)./p1(2:end)*100;
     
    hbar.iterate(1); % update progress by one iteration
end
time = toc;
close(hbar); % close the progress bar

%%%%%%%%%%%%%%%%%%%%%%    Convergence Criteria    %%%%%%%%%%%%%%%%%%%%%%%%%
% a) The RMSE should be less than tol.                                    %
% b) Variation in model parameters should be less one percent.            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
isConverged = rmse<1e-4 & all(err<1, 2);

if all(isConverged)
    fprintf('optimisation module works fine in class %s. (time per case = %1.1f milisecond)\n', modelName, time/N*1000);
else
    fprintf('optimissation converged for %d/%d tests for %s model. (time per case = %1.1f milisecond)\n', sum(isConverged), N, modelName, time/N*1000);
end

%% ellipticalCylinder
modelName = 'ellipticalCylinder';
nParams = 6;
nMultiRun = 10;

load(fullfile(path_to_data, strcat(modelName, '_rotated.mat')), 'testData', 'scheme');

nScheme = size(scheme, 1);
N    = length(testData);

flag = zeros(N, 1);
rmse = zeros(N, 1);
err  = zeros(N, nParams);

P1   = zeros(N, nParams);
P2   = zeros(N, nParams);

hbar = parfor_progressbar(N, ...
            sprintf('Optimisation is running for %s ...', modelName));
tic;        
parfor n = 1:N
    % read the input parameters
    p1 = [testData(n).diff, testData(n).r1, testData(n).r2,...
          testData(n).theta, testData(n).phi, testData(n).alpha]';
      
    P1(n,:) = p1;       
    
    % model fitting from random initial condition
    model = ellipticalCylinder('fitter',  optimizer('algorithm', algorithm));
    model.addConstraint('s0=1');
    
    %model.randomInit();
    %p = model.fit(scheme, testData(n).s);
    
    p = model.fitMultiRun(scheme, testData(n).s, nMultiRun);
    
    flag(n) = p(1);
    thisCost = p(2);
    rmse(n) = sqrt(thisCost/nScheme);
    
    p2 = p(4:end); 
    P2(n,:) = p2;
  
    err(n, :) = abs(p1-p2)./p1*100;
     
    hbar.iterate(1); % update progress by one iteration
end
time = toc;
close(hbar); % close the progress bar

%%%%%%%%%%%%%%%%%%%%%%    Convergence Criteria    %%%%%%%%%%%%%%%%%%%%%%%%%
% a) The RMSE should be less than tol.                                    %
% b) Variation in model parameters should be less one percent.            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
isConverged = rmse<1e-4 & all(err<1, 2);

if all(isConverged)
    fprintf('optimisation module works fine in class %s. (time per case = %1.1f milisecond)\n', modelName, time/N*1000);
else
    fprintf('optimissation converged for %d/%d tests for %s model. (time per case = %1.1f milisecond)\n', sum(isConverged), N, modelName, time/N*1000);
end

%% ball_stick
modelName = 'ball_stick';
nParams = 5;
nMultiRun = 10;

load(fullfile(path_to_data, strcat(modelName, '_rotated.mat')), 'testData', 'scheme');

nScheme = size(scheme, 1);
N    = length(testData);

flag = zeros(N, 1);
rmse = zeros(N, 1);
err  = zeros(N, nParams);

P1   = zeros(N, nParams);
P2   = zeros(N, nParams);

hbar = parfor_progressbar(N, ...
            sprintf('Optimisation is running for %s ...', modelName));
tic;        
parfor n = 1:N
    % read the input parameters
    s0 = 1;
    p1a =  testData(n).comp1_diff;
    p1b = [testData(n).comp2_diff, testData(n).comp2_theta, ...
           testData(n).comp2_phi];
    vFractions = [testData(n).comp1_vFraction, testData(n).comp2_vFraction]; 
    
    p1 = [vFractions(1), p1a, p1b]';
    P1(n,:) = p1;  
    
    comp1 = ball([s0, p1a]);
    comp2 = stick([s0, p1b]);
    
    % model fitting from random initial condition
    model = twoCompartmentModel(comp1, comp2, 'fitter',  optimizer('algorithm', algorithm));
    model.addConstraint('s0=1');
    
%     model.randomInit();
%     p = model.fit(scheme, testData(n).s);
    
    p = model.fitMultiRun(scheme, testData(n).s, nMultiRun);
    
    flag(n) = p(1);
    thisCost = p(2);
    rmse(n) = sqrt(thisCost/nScheme);
    
    p2 = p([4,7,9,10,11]); 
    P2(n,:) = p2;
  
    err(n, :) = abs(p1-p2)./p1*100;
     
    hbar.iterate(1); % update progress by one iteration
end
time = toc;
close(hbar); % close the progress bar

%%%%%%%%%%%%%%%%%%%%%%    Convergence Criteria    %%%%%%%%%%%%%%%%%%%%%%%%%
% a) The RMSE should be less than tol.                                    %
% b) Variation in model parameters should be less one percent.            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
isConverged = rmse<1e-4 & all(err<1, 2);

if all(isConverged)
    fprintf('optimisation module works fine in class %s. (time per case = %1.1f milisecond)\n', modelName, time/N*1000);
else
    fprintf('optimissation converged for %d/%d tests for %s model. (time per case = %1.1f milisecond)\n', sum(isConverged), N, modelName, time/N*1000);
end

%% cylinder_zeppelin
modelName = 'cylinder_zeppelin';
nParams = 9;
nMultiRun = 20;

load(fullfile(path_to_data, strcat(modelName, '_rotated.mat')), 'testData', 'scheme');

nScheme = size(scheme, 1);
N    = length(testData);

flag = zeros(N, 1);
rmse = zeros(N, 1);
err  = zeros(N, nParams);

P1   = zeros(N, nParams);
P2   = zeros(N, nParams);

hbar = parfor_progressbar(N, ...
            sprintf('Optimisation is running for %s ...', modelName));
tic;        
parfor n = 1:N
    % read the input parameters
    s0 = 1;
    p1a = [testData(n).comp1_diff , testData(n).comp1_r   , ...
           testData(n).comp1_theta, testData(n).comp1_phi];
    p1b = [testData(n).comp2_diffPar, testData(n).comp2_diffPerp, ... 
           testData(n).comp2_theta  , testData(n).comp2_phi];
    vFractions = [testData(n).comp1_vFraction, testData(n).comp2_vFraction]; 
    
    p1 = [vFractions(1), p1a, p1b]';
    P1(n,:) = p1;  
    
    comp1 = cylinder([s0, p1a]);
    comp2 = zeppelin([s0, p1b]);
    
    % model fitting from random initial condition
    model = twoCompartmentModel(comp1, comp2, 'fitter',  optimizer('algorithm', algorithm));
    model.addConstraint('s0=1');
    
    %model.randomInit();
    %p = model.fit(scheme, testData(n).s);
    
    p = model.fitMultiRun(scheme, testData(n).s, nMultiRun);
    
    flag(n) = p(1);
    thisCost = p(2);
    rmse(n) = sqrt(thisCost/nScheme);
    
    p2 = p([4, 7:10, 12:15]); 
    P2(n,:) = p2;
  
    err(n, :) = abs(p1-p2)./p1*100;
     
    hbar.iterate(1); % update progress by one iteration
end
time = toc;
close(hbar); % close the progress bar

%%%%%%%%%%%%%%%%%%%%%%    Convergence Criteria    %%%%%%%%%%%%%%%%%%%%%%%%%
% a) The RMSE should be less than tol.                                    %
% b) Variation in model parameters should be less one percent.            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
isConverged = rmse<1e-4 & all(err<1, 2);

if all(isConverged)
    fprintf('optimisation module works fine in class %s. (time per case = %1.1f milisecond)\n', modelName, time/N*1000);
else
    fprintf('optimissation converged for %d/%d tests for %s model. (time per case = %1.1f milisecond)\n', sum(isConverged), N, modelName, time/N*1000);
end

%% ellipticalCylinderGPD_ball_stick
modelName = 'ellipticalCylinder_ball_stick';
nParams = 12;
nMultiRun = 25;

load(fullfile(path_to_data, strcat(modelName, '_rotated.mat')), 'testData', 'scheme');

nScheme = size(scheme, 1);
N    = length(testData);

flag = zeros(N, 1);
rmse = zeros(N, 1);
err  = zeros(N, nParams);

P1   = zeros(N, nParams);
P2   = zeros(N, nParams);

hbar = parfor_progressbar(N, ...
            sprintf('Optimisation is running for %s ...', modelName));
tic;        
parfor n = 1:N
    % read the input parameters
    s0 = 1;
    p1a = [testData(n).comp1_diff, testData(n).comp1_r1   , ...
           testData(n).comp1_r2  , testData(n).comp1_theta, ...
           testData(n).comp1_phi , testData(n).comp1_alpha];
    p1b =  testData(n).comp2_diff;
    p1c = [testData(n).comp3_diff, testData(n).comp3_theta, ...
           testData(n).comp3_phi];
    vFractions = [testData(n).comp1_vFraction, ...
                  testData(n).comp2_vFraction, testData(n).comp3_vFraction]; 
    
    p1 = [vFractions(1:2), p1a, p1b, p1c]';
    P1(n,:) = p1;  
    
    comp1 = ellipticalCylinder([s0, p1a]);
    comp2 = ball([s0, p1b]);
    comp3 = stick([s0, p1c]);
    
    % model fitting from random initial condition
    model = threeCompartmentModel(comp1, comp2, comp3, 'fitter', optimizer('algorithm', algorithm));
    model.addConstraint('s0=1'); 
    data = testData(n).s;
    
%     model.randomInit();
%     p = model.fit(scheme, data);
    
    p = model.fitMultiRun(scheme, testData(n).s, nMultiRun);
    
    flag(n) = p(1);
    thisCost = p(2);
    rmse(n) = sqrt(thisCost/nScheme);
    
    p2 = p([4,5,8:13, 15, 17:19]); 
    P2(n,:) = p2;
  
    err(n, :) = abs(p1-p2)./p1*100;
     
    hbar.iterate(1); % update progress by one iteration
end
time = toc;
close(hbar); % close the progress bar

%%%%%%%%%%%%%%%%%%%%%%    Convergence Criteria    %%%%%%%%%%%%%%%%%%%%%%%%%
% a) The RMSE should be less than tol.                                    %
% b) Variation in model parameters should be less one percent.            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
isConverged = rmse<1e-4 & all(err<1, 2);

if all(isConverged)
    fprintf('optimisation module works fine in class %s. (time per case = %1.1f milisecond)\n', modelName, time/N*1000);
else
    fprintf('optimissation converged for %d/%d tests for %s model. (time per case = %1.1f milisecond)\n', sum(isConverged), N, modelName, time/N*1000);
end