function writeModel(obj, filename, varargin)
% saveModelConfig write the model configuration into a simple text file.

% Mohsen Farzi
% Email: m.farzi@leeds.ac.uk

% check if model is not empty
if isempty(obj.model)
    error('MATLAB:MYOSCOPE:writeModel',...
          'The model is not set yet.');
end

p = inputParser;
p.CaseSensitive = false;
p.addOptional('mode', 'new', @(v) assert(strcmp(v, 'new')||strcmp(v, 'append')));
p.parse(varargin{:});

mode = p.Results.mode;

[filepath,name,ext] = fileparts(filename);
if isempty(filepath)
    filepath = pwd;
end
if ~strcmp(ext, '.myo')
    warning('MATLAB:MYOSCOPE:writeParams',...
            'The filename extension is changed to ".myo".');
    ext = '.myo';    
end
filename = fullfile(filepath, strcat(name,ext));

if strcmp(mode, 'new')
    if isfile(filename)
        warning('MATLAB:MYOSCOPE:writeParams',...
                'Filename is overwritten.\n %s', filename);
    end
    fileId = fopen(filename,'w');
else
    fileId = fopen(filename,'a');
end

%% write model
fprintf(fileId, '##model\n');

% model name
fprintf(fileId, '#$name\n');
fprintf(fileId, '%s\n', obj.model.name);

fprintf(fileId, '#$paramsNum\n');
fprintf(fileId, '%d\n', obj.model.getParamsNum('all'));

fprintf(fileId, '#$freeParamsNum\n');
fprintf(fileId, '%d\n', obj.model.getParamsNum('free'));

fprintf(fileId, '\n');

% hyperparameters (name, value) pairs
nHyperparams = obj.model.getHyperparamsNum();
fprintf(fileId, '#$hyperparameters\n');
fprintf(fileId, '%d\n', nHyperparams);
if nHyperparams > 0
    hyperparams = obj.model.getHyperparams;
    hyperparamsList = obj.model.getHyperparamsList();
    for idx = 1:nHyperparams
        fprintf(fileId, '#$%s\n', hyperparamsList{idx});
        fprintf(fileId, '%d\n', hyperparams(idx));
    end
end
fprintf(fileId, '\n');

% constraints
constraintList = obj.model.getConstraints();
fprintf(fileId, '#$constraints\n');
fprintf(fileId, '%d\n', length(constraintList));
for idx = 1:length(constraintList)
    fprintf(fileId, '%s\n', constraintList{idx});
end
fprintf(fileId, '##end\n');
%\\
end % of writeModel