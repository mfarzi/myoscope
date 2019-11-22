function  property = searchConfig(fileName, varargin)
% SEARCHCONFIG search the model Config File 
%   property = SEARCHCONFIG(fileName, varargin) 
%   returns the properties of interest.
% 
%   flag returns an integer number for each property
%   0 --> Search is unsuccessful
%   1 --> A unique value is returned for the property
%   n --> n different matches exist for the property of interest

% Mohsen Farzi
% Email: m.farzi@leeds.ac.uk
%--------------------------------------------------------------------------
% check the number of arguments
if nargin<2
    error('MATLAB:loadCompartmentModel:InsufficientInputArgument', ...
          'Not enough input arguments.');
end
    
% check arguments
if ~isfile(fileName)
    error('MATLAB:loadCompartmentModel:fileIsNotFound', ...
          'File does not exist:\n%s', fileName);
end

% read config file
fileId = fopen(fileName, 'r');
result = textscan(fileId, '%s', 'delimiter', '\n');
modelConfig = result{1};
fclose(fileId);

% Allocate memory for output variables
nProps = nargin -1;
property = initProperty(varargin{:});

% search for each property
for idx = 1:nProps
    switch varargin{idx}
        case 'name'
            thisLineNo = cellfun(@(c) startsWith(c, 'name:'), modelConfig);
            if any(thisLineNo)
                tmp = textscan(modelConfig{thisLineNo}, 'name: %s');         
                modelName = tmp{1}{1};
                property.name = modelName;
            end
        %\\    
        case 'constraints'
            thisLine = cellfun(@(c) startsWith(c, 'constraints:'), modelConfig);
            thisLineNo = find(thisLine);
            tmp = textscan(modelConfig{thisLineNo}, 'constraints: %d');
            numberOfConstraints = tmp{1};
            property.constraints = modelConfig(thisLineNo+1:thisLineNo+numberOfConstraints);
        case 'params'    
            thisLine = cellfun(@(c) startsWith(c, 'parameters:'), modelConfig);
            thisLineNo = find(thisLine);
            property.params = sscanf(modelConfig{thisLineNo+2}, '%f');
        case 'hyperparams'
            thisLine = cellfun(@(c) startsWith(c, 'hyperparameters:'), modelConfig);
            thisLineNo = find(thisLine);
            tmp = textscan(modelConfig{thisLineNo}, 'hyperparameters: %d');
            numberOfHyperparams = tmp{1};
            property.hyperparams = modelConfig(thisLineNo+1:thisLineNo+numberOfHyperparams);
        %\\    
        otherwise
            error('MATLAB:searchConfig:unknownProperty', ...
                  'Property %s is not recognised.', varargin{idx});
    end
end
end

function property = initProperty(varargin)
% INITPROPERTY initialise the property files with empty matrix
fieldValuePair = cell(2*nargin, 1);
fieldValuePair(1:2:2*nargin) = varargin';
property = struct(fieldValuePair{:});
end