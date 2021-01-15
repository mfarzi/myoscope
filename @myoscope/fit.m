function obj = fit(obj, scheme, data, varargin) 
% Parse the parameter-value pairs
p = inputParser;
p.CaseSensitive = false;
p.addOptional('voxels', [], @(v) assert(isa(v, 'roi')));
p.addParameter('nMultiRun'   , 20, @(v) assert(isscalar(v) && v>0 && rem(v,1)==0));


p.parse(varargin{:});              

% update the fitter
nMultiRun = p.Results.nMultiRun;
voxels = p.Results.voxels;

nPxls = size(data, 2);
params = zeros(obj.model.getParamsNum+2, nPxls);

% to measure parfor loop pgrogress and temperory storage 
tmpFolder = fullfile(obj.cwd, 'tmp');
if ~isfolder(tmpFolder)
    mkdir(tmpFolder);
else
    delete(fullfile(tmpFolder, '*.myo'));
end

thisModel = obj.model;
fprintf('start myoscope %s \n', thisModel.name);
tic;
parfor i = 1:nPxls
    p = thisModel.fitMultiRun(scheme, data(:,i), nMultiRun);
    params(:,i) = p;
    id = fopen(fullfile(tmpFolder, sprintf('%d.myo', i)), 'w');
    fprintf(id, '%1.6e\n', p);
    fclose(id);
end
runTime = toc;
fprintf('end myoscope-run time=%f minutes \n', runTime/60);
obj.params = params(3:end,:);

% compute RMSE
nScheme = size(scheme,1);
squaredError = params(2,:);
rmse = sqrt(squaredError/nScheme);
obj.rmse = rmse;

% compute AIC
sigma = rmse;
logLikelihood = -squaredError./(2*sigma.^2) - nScheme*log(sqrt(2*pi)*sigma);
aic = 2*thisModel.getParamsNum('free') - 2*logLikelihood;
obj.aic = aic;

obj.runTime = runTime;
obj.voxels = voxels.