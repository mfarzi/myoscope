function fileName = fit_zeppelin_tensor_ball_02(paramsFolder, ...
                    path2voxelData, forceWrite)
% FIT_ZEPPELIN_TENSOR_BALL_01 constitues the forth step in the biophysical 
% modeling pipeline 
%
% DESCRIPTION: fits a prolate, an oblate, and a spherical ellipsoids to
% data to explain anistropy as decomposition of isotropic extracellular
% space, intracellular with parralel myocytes, and intracellular space with
% myocytes uniformly distributed in the Azimuth angle.
%
% Author: Mohsen Farzi                                                    
% Email:  m.farzi@leeds.ac.uk                                             
%\\

[~, dataName, ~] = fileparts(path2voxelData);
thisName = strrep(dataName, 'data', 'params');
thisName = strcat(thisName, '.mat');

modelName = 'zeppelin-tensor-ball-02';
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
nParams = 20;
nMultiRun = 50;

%% fit model
params = zeros(nPxlTot, nParams);
hbar = parfor_progressbar(nPxlTot, 'Optimisation is running ...');
scheme = scheme; % make it accessible on the workers when using parfor
tic;
parfor thisPxl = 1:nPxlTot
    data = s(:,thisPxl);

    cmp1 = zeppelin();                      % prolate ellipsoid
    cmp2 = tensor();                        % oblate spheroid
    cmp3 = ball();                          % ball

    model = threeCompartmentModel(cmp1, cmp2, cmp3);
    model.addConstraint('s0 = 1');
    
    model.addConstraint('comp2_theta = comp1_theta');
    model.addConstraint('comp2_phi   = comp1_phi'  );

    model.addConstraint('comp3_diff >= 1.8e-9');
       
    model.addConstraint('comp2_diffPar=(comp1_diffPar+comp1_diffPerp)*0.5');
    model.addConstraint('comp2_diffPerp1=(comp1_diffPar+comp1_diffPerp)*0.5');
    model.addConstraint('comp1_diffPerp >= comp2_diffPerp2');
    
%     model.randomInit();
%     p = model.fit(scheme, data);
    
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