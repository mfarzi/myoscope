function model = load(fileName)
% LOADCOMPARTMENTMODEL load compartment model 
%   model = LOADCOMPARTMENTMODEL(filename)
%   returns a multi-compartment object with given configurations.
%
% see also compartment

% Mohsen Farzi
% Email: m.farzi@leeds.ac.uk

% check arguments
if ~isfile(fileName)
    error('MATLAB:loadCompartmentModel:fileIsNotFound', ...
          'File does not exist:\n%s', fileName);
end

% read config file
config = searchConfig(fileName, 'name', 'constraints', 'params', 'hyperparams');

% construct the model
model = str2model(config.name);

% make sure links are of the same type
tmp = textscan(config.constraints{1}, '%s', 'delimiter', {':',','});
linkType = tmp{1};
linkType(1) = [];

numOfLinks = length(model.links);
for thisLinkNo = 1:numOfLinks
    isDummy = strcmp(linkType{thisLinkNo}, 'dummy');
    isMatched = strcmp(linkType{thisLinkNo}, model.links(thisLinkNo).type);
    if ~isMatched && ~isDummy
        model.links(thisLinkNo).set('type', linkType{thisLinkNo});
    end
end

% check for other constraints
defaultConstraints = model.getConstraints;

idx = ismember(config.constraints, defaultConstraints);
idx(1) = true; % link types have been checked previously

thisConstraints = config.constraints(~idx);
% enforce constraints
nConstraints = length(thisConstraints);
for i=1:nConstraints
    model.addConstraint(thisConstraints{i});
end

% set hyperparameters
nHyperparams = length(config.hyperparams);
hparams = zeros(nHyperparams, 1);
for i = 1:nHyperparams
    tmp = textscan(config.hyperparams{i},'%s %f');
    hparams(i) = tmp{2};
end
if ~isempty(hparams)
    model.set('hyperparams', hparams);
end

% set parameters values
model.set('params', config.params);
end % of loadCompartmentModel


    

