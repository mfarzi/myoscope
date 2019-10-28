function myoscopeSaveModel(model, fileName)
% saveModelConfig write the model configuration into a simple text file.

% Mohsen Farzi
% Email: m.farzi@leeds.ac.uk

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
fprintf(fileId, 'hyperparameters: %d\n', model.getHyperparamsNum());
nHyperparams = model.getHyperparamsNum();
if nHyperparams > 0
    hyperparams = model.getHyperparams;
    hyperparamsList = model.getHyperparamsList();
    for idx = 1:nHyperparams
        fprintf(fileId, '%s: %d\n', hyperparamsList{idx}, hyperparams(idx));
    end
end
fprintf(fileId, '\n');

% constraints
constraintList = model.getConstraints();
fprintf(fileId, 'constraints: %d\n', length(constraintList));
for idx = 1:length(constraintList)
    fprintf(fileId, '%s\n', constraintList{idx});
end
fprintf(fileId, '\n');

% parameters
fprintf(fileId, 'parameters: %d\n', model.getParamsNum());
params = model.getParams;
paramsList = model.getParamsList();
nParams = model.getParamsNum();
% print variable names in the first row
for idx = 1:nParams
    fprintf(fileId, '%s ', paramsList{idx});
end
fprintf(fileId, '\n');
% print parameter values in the second row
for idx = 1:nParams
    fprintf(fileId, '%1.6e ', params(idx));
end
%\\
end % of saveModelConfig