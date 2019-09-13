% This script test if the optimisation procedure converges to the given
% parameters for each model.

clc; clear all; close all
addpath('../models');
addpath('../');
addpath('../linker');
addpath('../utility');
addpath('../../utility/');

path_to_data = ['./testData/synthetic/'];

algorithm = 'levenberg-marquardt'; %'interior-point' %'conjugate-gradient';            
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
    model.addConstraint('diffPar>=diffPerp1');
    model.addConstraint('diffPerp1>=diffPerp2');
    
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

%% ball with Tensor
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
    model = tensor('fitter', optimizer('algorithm', algorithm));
    model.addConstraint('s0=1');
    model.addConstraint('diffPar=diffPerp1');
    model.addConstraint('diffPar=diffPerp2');
    
    %model.randomInit();
    %p = model.fit(scheme, testData(n).s);
    
    p = model.fitMultiRun(scheme, testData(n).s, nMultiRun);
    
    flag(n) = p(1);
    thisCost = p(2);
    rmse(n) = sqrt(thisCost/nScheme);
    
    p2 = p(4); 
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

%% zeppelin with Tensor
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
    model = tensor('fitter', optimizer('algorithm', algorithm));
    model.addConstraint('s0 = 1');
    model.addConstraint('diffPerp1 = diffPerp2');
    model.addConstraint('diffPar >= diffPerp1');
    
    %model.randomInit();
    %p = model.fit(scheme, testData(n).s);
    
    p = model.fitMultiRun(scheme, testData(n).s, nMultiRun);
    
    flag(n) = p(1);
    thisCost = p(2);
    rmse(n) = sqrt(thisCost/nScheme);
    
    if p(4) == p(5)
        p2 = p([5,6,7,8]);
    else
        p2 = p([4,5,7,8]);
    end
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

%% stick with tensor
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
    model = tensor('fitter',  optimizer('algorithm', algorithm));
    model.addConstraint('s0=1');
    model.addConstraint('diffPerp1=0');
    model.addConstraint('diffPerp2=0');
    
    %model.randomInit();
    %p = model.fit(scheme, testData(n).s);
    
    p = model.fitMultiRun(scheme, testData(n).s, nMultiRun);
    
    flag(n) = p(1);
    thisCost = p(2);
    rmse(n) = sqrt(thisCost/nScheme);
    
    p2 = p([4,7,8]); 
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

%% ball_stick with tensor
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
    
    comp1 = tensor('fitter', optimizer('algorithm', algorithm));
    comp1.addConstraint('s0=1');
    comp1.addConstraint('diffPar=diffPerp1');
    comp1.addConstraint('diffPar=diffPerp2');
    comp1.addConstraint('theta=0');
    comp1.addConstraint('phi=0');
    comp1.addConstraint('alpha=0');
    
    % model fitting from random initial condition
    comp2 = tensor('fitter',  optimizer('algorithm', algorithm));
    comp2.addConstraint('s0=1');
    comp2.addConstraint('diffPerp1=0');
    comp2.addConstraint('diffPerp2=0');
    comp2.addConstraint('alpha=0');
    
    
    % model fitting from random initial condition
    model = twoCompartmentModel(comp1, comp2, 'vFractions', vFractions, 's0', 1, 'fitter',  optimizer('algorithm', algorithm));
    model.addConstraint('s0=1');
    
%     model.randomInit();
%     p = model.fit(scheme, testData(n).s);
    
    p = model.fitMultiRun(scheme, testData(n).s, nMultiRun);
    
    flag(n) = p(1);
    thisCost = p(2);
    rmse(n) = sqrt(thisCost/nScheme);
    
    p2 = p([4, 7, 14, 17, 18]); 
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