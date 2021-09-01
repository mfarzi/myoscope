function writeModel(obj, filename, varargin)
% writeModel is a method for class MULTICOMPARTMENT
% 
% writeModel(obj, filename) save the model configuration into a text file.
%             Input args:
%               filename: Filename of type char to write the result.
%
% writeModel(obj, filename, permission) save the model configuration with 
% the type of access specified by permission.
%             Input args:
%               filename: Filename of type char to write the result.
%             permission: 'w' (default) | 'a'
%                         'w': creates a new file. Discard existing
%                              contents, if any.
%                         'a': opens or creates a new file. Append data to 
%                              the end of the file.
%
% writeModel(...,<param>,<value>,...) allows passing more information as
% param-value pairs.
%             parameters:
%              overwrite: Logical flag to whether write over existing file 
%                         [Default=false]               
%
% Mohsen Farzi
% Email: m.farzi@leeds.ac.uk

% check input arguments
filename = validateFilename(filename);
permission = parseInputVariables(varargin{:});
if strcmp(permission, 'w')
    fileId = fopen(filename,'w');
else
    fileId = fopen(filename,'a');
end

%% write model
fprintf(fileId, '##model\n');

% model name
fprintf(fileId, '#$name 1 string\n');
fprintf(fileId, '%s\n', obj.name);

% parameter names
fprintf(fileId, '#$paramsName 1 string\n');
fprintf(fileId, '%s ', obj.getParamsName{:});
fprintf(fileId, '\n');

% hyperparameters (name, value) pairs
if obj.nHyperparams > 0
    [hparams, hparamsList] = obj.getHyperparams();
    fprintf(fileId, '#$hyperparamsName 1 string\n');
    fprintf(fileId, '%s ', hparamsList{:});
    fprintf(fileId, '\n');
    
    fprintf(fileId, '#$hyperparams 1 numeric\n');
    fprintf(fileId, '%d ', hparams);
    fprintf(fileId, '\n');
end

% constraints
constraintList = obj.getConstraints();
nConstraints = length(constraintList);
fprintf(fileId, sprintf('#$constraints %d string\n', nConstraints));
for n = 1:nConstraints
    fprintf(fileId, '%s\n', constraintList{n});
end
fprintf(fileId, '##ENDmodel\n');
fclose(fileId);
%\\
end % of writeModel

function permission = parseInputVariables(varargin)
    p = inputParser;
    p.CaseSensitive = false;
    p.addOptional('permission', 'w', ...
        @(v) assert(isa(v, 'char') && (strcmp(v, 'a')||strcmp(v, 'w')),...
        'MATLAB:MULTICOMPARTMENT:writeModel',...
        "Possible permissions include 'a' or 'w'."));
    p.parse(varargin{:});

    permission = p.Results.permission;
end