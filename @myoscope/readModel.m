function obj = readModel(filename)
% LOADCOMPARTMENTMODEL load compartment model 
%   model = LOADCOMPARTMENTMODEL(filename)
%   returns a multi-compartment object with given configurations.
%
% see also compartment

% Mohsen Farzi
% Email: m.farzi@leeds.ac.uk

% check arguments
[filepath,name,ext] = fileparts(filename);
if isempty(filepath)
    filepath = pwd;
    filename = fullfile(filepath, strcat(name, ext));
end
if ~isfile(filename)
    error('MATLAB:MYOSCOPE:readModel', ...
          'File does not exist:\n%s', filename);
end

%%
modelExist = false;
fileId = fopen(filename, 'r');
thisLine = fgetl(fileId);
while ischar(thisLine)
    if startsWith(thisLine, '##model')
        modelExist = true;
        break;
    end
    thisLine = fgetl(fileId);
end

if not(modelExist)
    error('MATLAB:MYOSCOPE:readModel',...
          'Header information "##model" cannot be found.');
end
%% read model config
config = struct;
while ischar(thisLine) || strcmp(thisLine, '##end')
    switch thisLine
        case '#$name'
            fieldVal = fgetl(fileId);
            config.name = fieldVal;
            
        case '#$constraints'
            constraintsNo = str2double(fgetl(fileId));
            fieldVal = cell(constraintsNo,1);
            for i=1:constraintsNo
                fieldVal{i} = fgetl(fileId);
            end
            config.constrains = fieldVal;
            
        case '#$hyperparameters'
            config.hyperparamsNum = str2double(fgetl(fileId));
            if config.hyperparamsNum>0
                hyperparams = zeros(hyperparamsNum, 1);
                hyperparamsList = cell(hyperparamsNum, 1);
                for i=1:hyperparamsNum
                    hyperparamsList{i} = fgetl(fileId);
                    hyperparams(i) = str2double(fgetl(fileId));
                end
                config.hyperparams = hyperparams;
                config.hyperparamsList = hyperparamsList;
            end
    end%switch-case
    
    thisLine = fgetl(fileId);
end
fclose(fileId);
            
%% construct the model
obj = myoscope(config.name);

% make sure links are of the same type
tmp = textscan(config.constraints{1}, '%s', 'delimiter', {':',','});
linkType = tmp{1};
linkType(1) = [];

numOfLinks = length(obj.model.links);
for thisLinkNo = 1:numOfLinks
    isDummy = strcmp(linkType{thisLinkNo}, 'dummy');
    isMatched = strcmp(linkType{thisLinkNo}, obj.model.links(thisLinkNo).type);
    if ~isMatched && ~isDummy
        obj.model.links(thisLinkNo).set('type', linkType{thisLinkNo});
    end
end

% check for other constraints
defaultConstraints = obj.model.getConstraints;

idx = ismember(config.constraints, defaultConstraints);
idx(1) = true; % link types have been checked previously

thisConstraints = config.constraints(~idx);
% enforce constraints
nConstraints = length(thisConstraints);
for i=1:nConstraints
    obj.model.addConstraint(thisConstraints{i});
end

% set hyperparameters
if config.hyperparamsNum>0
    obj.model.set('hyperparams', config.hyperparams);
end
end%of readModel


    

