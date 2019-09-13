function fileName = fit_cylinderGPD_ball_02(paramsFolder, path2voxelData, ...
                    forceWrite)
% FIT_CYLINDER_BALL_01 constitues the forth step in the biophysical 
% modeling pipeline 
%
% Author: Mohsen Farzi                                                    
% Email:  m.farzi@leeds.ac.uk                                             
%\\

[~, dataName, ~] = fileparts(path2voxelData);
thisName = strrep(dataName, 'data', 'params');
thisName = strcat(thisName, '.mat');

modelName = 'cylinderGPD-ball-02';
thisFolder = fullfile(paramsFolder, modelName);
if ~isfolder(thisFolder)
    mkdir(thisFolder);
end

fileName = fullfile(thisFolder, thisName);

if isfile(fileName) && not(forceWrite)
    fprintf(strcat('Skipped step 4. Params exists for %s model.', ...
                   '\nSee %s \n', ...
                   'If you wish to save a different copy,',...
                   ' change "expNo".\nIf you wish to overwrite', ...
                   ' the existing file, turn the flag "forceWrite"',...
                   ' true.\n'), modelName, fileName);
   return;
else
    fprintf('Step 4: fittimg params for %s model.\n', modelName);
end

%% load data
load(path2voxelData, 's', 'scheme', 'mask');
nPxlTot = size(s, 2);
nParams = 12;
nMultiRun = 50;

%% fit model
params = zeros(nPxlTot, nParams);
hbar = parfor_progressbar(nPxlTot, 'Optimisation is running ...');
scheme = scheme; % make it accessible on the workers when using parfor
tic;
parfor thisPxl = 1:nPxlTot
    data = s(:,thisPxl);

    cmp1 = cylinderGPD(); % cylinder with Guassian Phase Distribution
    cmp2 = ball();        % ball

    model = twoCompartmentModel(cmp1, cmp2);
    model.addConstraint('s0 = 1');
    model.addConstraint('comp2_diff >= 1.8e-9');
    
    
    % model.randomInit();
    % p = model.fit(Scheme, data);
    
    p = model.fitMultiRun(scheme, data, nMultiRun);
    
    params(thisPxl, :) = p;
    
    fprintf(strcat('pixel %d is done with cost %d and exit flag %d', ...
                   'and muscle fraction %1.2f.\n'), thisPxl, p(2), ...
                   p(1), p(4)*100);
               
    hbar.iterate(1); % update progress by one iteration
end
close(hbar); % close the progress bar
runTime = toc;

% save results
save(fileName, 'params', 'mask', 'runTime', 'modelName');