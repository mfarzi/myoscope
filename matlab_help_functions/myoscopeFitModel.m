function myoscopeFitModel(path2model, path2data, path2scheme, path2output, varargin)
% MYOSCOPEFITMODEL fit compartment model to DW-MRI data 
%   MYOSCOPEFITMODEL use numerical optimisation techniques to fit model
%   parameters to the given diffusion weighted MR data.
%
% Input:
%           path2model: Compartmental Model
%
%            path2data: Reconstructed dw-mri data
%
%          path2scheme: Stejeskal-Tanner diffusion scheme 
%
%          path2output: Output folder to store fitted parameters
%
%           Name,Value: [optional] Using name-value pairs to control
%                       different aspects of the function
%
%                           "prettyPrint": Print Progress for each voxel
%                                          Default: false
%
%                              "fileName": output file name.
%                                          Default: parameters
%
%                                   "roi": path to a ROI data
%                                          Default: ''
%
% see also compartment

% Mohsen Farzi
% Email: m.farzi@leeds.ac.uk

% check input variables
p = inputParser;
p.CaseSensitive = false;
p.addParameter('prettyPrint'   , false, ...
    @(v) assert(islogical(v)&&isscalar(v), 'Value must be logical.'));
p.addParameter('fileName'   , 'parameters', ...
    @(v) assert(ischar(v), 'Value must be of type char.'));
p.addParameter('roi'   , '', ...
    @(v) assert(ischar(v), 'Value must be of type char.'));
%\\
p.parse(varargin{:});              
prettyPrint = p.Results.prettyPrint;
fileName = p.Results.fileName;
path2roi = p.Results.roi;

if ~isempty(path2roi) && ~isfile(path2roi)
    error('MATLAB:myoscopeFitModel:fileIsNotFound', ...
          'File does not exist:\n%s', path2roi);
end

if ~isfile(path2model)
    error('MATLAB:myoscopeFitModel:fileIsNotFound', ...
          'File does not exist:\n%s', path2model);
end

if ~isfile(path2data)
    error('MATLAB:myoscopeFitModel:fileIsNotFound', ...
          'File does not exist:\n%s', path2data);
end

if ~isfile(path2scheme)
    error('MATLAB:myoscopeFitModel:fileIsNotFound', ...
          'File does not exist:\n%s', path2scheme);
end

if ~isfolder(path2output)
    error('MATLAB:myoscopeFitModel:folderIsNotFound', ...
          'Output folder does not exist:\n%s', path2output);
end

% read model
model = myoscopeLoadModel(path2model);

% read diffusion data
data = myoscopeLoadData(path2data);

% read scheme file
scheme = myoscopeReadScheme(path2scheme);
% ad hoc solution for gradient directions
for i = 1:61:305
    scheme.G_dir(i) = 0;
    scheme.G_dir(i+1:i+60) = [1:10, 1:10, 1:10, 1:10, 1:10, 1:10];
end

% initialize params matrix
nPxl = size(data, 2);
params = zeros(model.getParamsNum + 2, nPxl);
% fit model to data
parfor thisPxl = 1:5%nPxl
    p = model.fitMultiRun(scheme, data(:,thisPxl));
    params(:,thisPxl) = p;
    if prettyPrint
        fprintf('pixel %d/%d: exit flag = %d, and cost function = %d \n',...
                thisPxl, nPxl, p(1), p(2));
    end
end

% write output results
path2result = fullfile(path2output, strcat(fileName, '.csv'));
fileId = fopen(path2result, 'w');
varNameList = [{'flag', 'cost'}, model.getParamsList()];
nParams = length(varNameList);
formatSpec = [repmat('%s,', 1, nParams-1), '%s\n'];
fprintf(fileId, formatSpec, varNameList{:});
formatSpec = [repmat('%1.6e,', 1, nParams-1), '%1.6e\n'];
fprintf(fileId, formatSpec, params);
fclose(fileId);
%\\
end % of myoscopeFitModel