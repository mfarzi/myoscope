function saveModelConfig(model, path2folder, name)
% saveModelConfig write the model configuration into a simple text file.

% Mohsen Farzi
% Email: m.farzi@leeds.ac.uk

fileName = fullfile(path2folder, name);
fileId = fopen(fileName, 'w');

% file header
fprintf(fileId, '# Multi-Compartment Model Structure Type\n');
fprintf(fileId, '# Version: v00.1\n');
fprintf(fileId, '\n');

% model name
fprintf(fileId, 'name: %s\n', model.name);
fprintf(fileId, 'number of parameters: %d\n', model.getParamsNum('all'));
fprintf(fileId, 'number of free parameters: %d\n', ...
                 model.getParamsNum('free'));
fprintf(fileId, 'number of hyperparameters: %d\n', ...
                 model.getHyperparamsNum());
fprintf(fileId, '\n');

% hyperparameters (name, value) pairs
fprintf(fileId, 'hyperparameters:\n');
nHyperparams = model.getHyperparamsNum();
if nHyperparams > 0
    hyperparams = model.getHyperparams;
    hyperparamsList = model.getHyperparamsList();
    for idx = 1:nHyperparams
        fprintf(fileId, '%s: %d\n', hyperparamsList{idx}, hyperparams(idx));
    end
else
    fprintf(fileId, 'NA');
end
fprintf(fileId, '\n');

% constraints
fprintf(fileId, 'constraints:\n');
constraintList = model.getConstraint();
for idx = 1:length(constraintList)
    fprintf(fileId, '%s\n', constraintList{idx});
end
fprintf(fileId, '\n');

% parameters
fprintf(fileId, 'parameters:\n');
params = model.getParams;
paramsList = model.getParamsList();
nParams = model.getParamsNum();
for idx = 1:nParams
    fprintf(fileId, '%s: %1.6e\n', paramsList{idx}, params(idx));
end

end             