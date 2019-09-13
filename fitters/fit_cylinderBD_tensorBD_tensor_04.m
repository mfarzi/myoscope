function fileName = fit_cylinderBD_tensorBD_tensor_04(paramsFolder, ...
                    path2voxelData, forceWrite)
% FIT_CYLINDERGDR_BALL_01 constitues the forth step in the biophysical 
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

modelName = '[[cylinderBD-tensorBD]-tensor]';
folderName = strcat(modelName,'-01');
tmpModel = str2model(modelName);

thisFolder = fullfile(paramsFolder, folderName);
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
nParams = tmpModel.getParamsNum() + 2;         % exist flag + cost function 
nMultiRun = 20;

%% fit model
params = zeros(nPxlTot, nParams);
hbar = parfor_progressbar(nPxlTot, 'Optimisation is running ...');
scheme = scheme; % make it accessible on the workers when using parfor
tic;
for thisPxl = 123%:nPxlTot
%for thisPxl = 148
    data = s(:,thisPxl);

    model = str2model(modelName);
    model.comp{1}.comp{1}.set('nBinsTheta', 50);
    model.comp{1}.comp{1}.set('nBinsPhi', 50);
    model.comp{1}.comp{2}.set('nBinsTheta', 50);
    model.comp{1}.comp{2}.set('nBinsPhi', 50);
    
    model.addConstraint('s0 = 1');
    
    model.addConstraint('comp1.comp1.kappa2 >= comp1.comp1.kappa1');
    model.addConstraint('comp1.comp2.kappa2 = comp1.comp1.kappa2');
    model.addConstraint('comp1.comp2.kappa1 = comp1.comp1.kappa1');
    
    model.addConstraint('comp1.comp2.phi = comp1.comp1.phi');
    model.addConstraint('comp1.comp2.theta = comp1.comp1.theta');
    model.addConstraint('comp1.comp2.alpha = comp1.comp1.alpha');
    
    model.addConstraint('comp1.comp2.diffPar = comp2.diffPar');
    model.addConstraint('comp1.comp2.diffPerp= comp2.diffPar*(1-comp1.comp1.s0)');
    
    model.addConstraint('comp2.phi = comp1.comp1.phi');
    model.addConstraint('comp2.theta = comp1.comp1.theta');
    model.addConstraint('comp2.alpha = comp1.comp1.alpha');
    model.addConstraint('comp2.diffPar >= 1.8e-9');
    model.addConstraint('comp2.diffPar >= comp1.comp1.diffPar');
    model.addConstraint('comp2.diffPerp1 = comp2.diffPar');
    model.addConstraint('comp2.diffPar >= comp2.diffPerp2');
    
%     model.addConstraint('comp1.comp1.s0>=0.8');
%     model.randomInit();
%     p = model.fit(scheme, data);
%     
    p = model.fitMultiRun(scheme, data, nMultiRun);
%     
    params(thisPxl, :) = p;
    
    fprintf(strcat('pixel %d is done with cost %d and exit flag %d', ...
                   'and sheetlet fraction %1.2f.\n'), thisPxl, p(2), ...
                   p(1), p(4)*100);
               
    hbar.iterate(1); % update progress by one iteration
end
close(hbar); % close the progress bar
runTime = toc;

% save results
save(fileName, 'params', 'mask', 'runTime', 'modelName');