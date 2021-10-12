function model = readModel(filename)
% readModel is a method for class MULTICOMPARTMENT
% 
% readModel(filename) return a MULTICOMPARTMENT object with input model
% configuration stored in filename.
%             Input args:
%               filename: Filename of type char to read the model
%            Output args:
%                  model: A multicompartment object
%
% Mohsen Farzi
% Email: m.farzi@leeds.ac.uk

% read model configurations
config = readConfig(filename, 'model');

%% construct the model
model = compartment.str2model(config.name{1}{1});
model.setParamsName(config.paramsName{1});

% remove default constraints which are not in the config
defaultConstraints = model.getConstraints;
configConstraints = cellfun(@(c) c{1}, config.constraints,...
    'UniformOutput', false);
nDeafultConstraints = length(defaultConstraints);
for i=1:nDeafultConstraints
    if not(ismember(defaultConstraints{i},configConstraints))
        % find the left-hand variable
        str = defaultConstraints{i};
        lhv = cellfun(@(c) startsWith(str, c), config.paramsName{1});
        % remove the active constraint
        if ~strcmp(model.links(lhv).type, 'independent')
            model.removeConstraint(config.paramsName{1}{lhv});
        end
    end
end

% enforce new constraints [not set by deafult]
nConstraints = length(config.constraints);
defaultConstraints = model.getConstraints;
for i=1:nConstraints
    if not(ismember(config.constraints{i}{1}, defaultConstraints))
        model.addConstraint(config.constraints{i}{1});
    end
end

% set hyperparameters
if model.nHyperparams>0
    model.hyperparams = config.hyperparams{1};
    model.hyperparamsName = config.hyperparamsName{1};
end
end%of readModel
