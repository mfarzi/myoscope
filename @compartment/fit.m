function [params, rmse, exitFlag] = fit(obj, data, schemefile, varargin)
% fit is a method for class MULTICOMPARTMENT
%
% fit(obj, data, scheme) fit parameter models to the obererved data.
%             Input args:
%                   data: A double matrix of MxN. All values must be 
%                         normalised to the b0-signal. M is the number DW
%                         measurements and N is the number of voxels.
%             schemefile: A scheme object with M measurements.
%            Output args:
%                      p: A double matrix of size PxN. P equals the number
%                         of model parameters. 
%                   rmse: Root Mean Squared Error
%               exitFlag: Return the exit flag for the optimisation
%
% fit(...,<param>,<value>,...) allows passing more information as
% param-value pairs.
%             parameters:
%                  nreps: Number of different runs with a random initial 
%                         point per voxel. [Default=1]
%                tfolder: Temporary folder to store intermediate results               
%

% check consistency between data and scheme file
assert(size(data, 1)==schemefile.measurementsNum(),...
    'MATLAB:multicompartment:fit',...
    'Number of rows must be equal in data and scheme.')

% parse input parameters
[nreps, tfolder, nWorkers, progress] = parseInputVariables(varargin{:});
nScheme = schemefile.measurementsNum();%size(schemefile, 1);
nVoxels = size(data, 2);

% setup the optimizer handle functions
hparams = obj.getHyperparams();
fitter = optimizer(); % an "optimizer" object for fitting model parameters
fitter.setfeval(@(x, s) obj.Feval(x, s, schemefile))
fitter.setfjac(@(x,~) obj.Fjac(x, schemefile));

% extend data by repeating each column 'nreps' times
N = nVoxels*nreps;
repIds = reshape(repmat(1:nVoxels, nreps, 1),1,N);
data = data(:,repIds);

% allocate memory to store results
params = zeros(obj.nParams, N);
rmse = zeros(1, N);
exitFlag = zeros(1,N);

% create working directory to store intermediate files
wrkdir = fullfile(tfolder, datestr(datetime,'yyyymmdd-HHMM'));
if progress
    mkdir(wrkdir);
end

parfor(n=1:N, nWorkers)
    % input signal
    sig = data(:, n);
    % generate initial values for the model parameters
    x0 = obj.links.initx();
    % fit model to sig
    [x, cost, flag] = fitter.run(x0, sig);
    p = obj.links.map(x);
    err = sqrt(cost/nScheme);
    % store results
    params(:,n) = p;
    rmse(n) = err;
    exitFlag(n) = flag;
    % print results
    if progress
        fileId = fopen(fullfile(wrkdir, sprintf('%d.myo', n)),'w');
        fprintf(fileId, '%1.6e\n', [flag; err; p]);
        fclose(fileId);
    end
end%of parfor

% select the min rmse per repetitions
if nreps>1
    rmse = reshape(rmse, [nreps, nVoxels]);
    [rmse, idx] = min(rmse);
    idx = idx+(0:nVoxels-1)*nreps;
    params = params(:,idx);
    exitFlag = exitFlag(idx);
end
end

function [nreps, tfolder, nWorkers, progress] = parseInputVariables(varargin)
% Parse the parameter-value pairs
p = inputParser;
p.CaseSensitive = false;
p.addParameter('nreps', 1, ...
    @(v) assert(mod(v,1)==0 && isscalar(v) && isnumeric(v) && v>0,...
    'MATLAB:multicompartment:fit',...
    'nreps must be a scalar integer value.'));
p.addParameter('tfolder'   , pwd,...
    @(v) assert(ischar(v) && isfolder(v),...
    'MATLAB:multicompartment:fit',...
    'tfolder must be a valid path to a directory.'));

p.addParameter('runInSerial'   , false,...
    @(v) assert(islogical(v) && isscalar(v),...
    'MATLAB:multicompartment:fit',...
    'runInsSerial must be a scalar logical variable.'));

p.addParameter('nWorkers', 12, ...
    @(v) assert(mod(v,1)==0 && isscalar(v) && isnumeric(v) && v>=0,...
    'MATLAB:multicompartment:fit',...
    'nreps must be a scalar, non-negative integer value.'));

p.addParameter('progress'   , false,...
    @(v) assert(islogical(v) && isscalar(v),...
    'MATLAB:multicompartment:fit',...
    'progress must be a scalar logical variable.'));

p.parse(varargin{:});              

% update the fitter
nreps = p.Results.nreps;
tfolder = p.Results.tfolder;
if p.Results.runInSerial
    nWorkers = 0;
else
    nWorkers = p.Results.nWorkers;
end
progress = p.Results.progress;
end