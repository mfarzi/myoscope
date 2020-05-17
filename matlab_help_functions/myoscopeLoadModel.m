function model = myoscopeLoadModel(fileName)
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
linkType = textscan(config.constraints{1}, '%s', 'delimiter', {':',','});
linkType{1} = [];

numOfLinks = length(model.links);
for thisLinkNo = 1:numOfLinks
    if strm(linkType{thisLinkNo}
    model.links(thisLinkNo).set('type', linkType{thisLinkNo});
end
defaultConstraints = model.getConstraints;

idx = ismember(config.constraints, defaultConstraints);
thisConstraints = config.constraints(~idx);
% enforce constraints
nConstraints = length(thisConstraints);
for i=1:nConstraints
    model.addConstraint(thisConstraints{i});
end

% set hyperparameters
%nHyperparams = length(config.hyperparams);
%for i = 1:nHyperparams
%    model.set(config.hyperparams{i});
%end

% set parameters values
model.set('params', config.params);
end % of loadCompartmentModel


    
