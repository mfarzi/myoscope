function  [property, flag] = searchConfig(fileName, varargin)
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
flag = zeros(nProps, 1);
property = initProperty(varargin{:});

% search for each property
for idx = 1:nProps
    switch varargin{idx}
        case 'name'
            thisLineNo = cellfun(@(c) startsWith(c, 'name:'), modelConfig);
            if any(thisLineNo)
                result = textscan(modelConfig{thisLineNo}, 'name: %s');         
                flag(idx) = length(result{1});
                property = setfield(property, 'name', result{1});
            end
            
        otherwise
            error('MATLAB:searchConfig:unknownProperty', ...
                  'Property %s is not recognised.', propName);
    end
end
end

function property = initProperty(varargin)
% INITPROPERTY initialise the property files with empty matrix
fieldValuePair = cell(2*nargin, 1);
fieldValuePair{1:2:2*nargin} = varargin{:};
property = struct(fieldValuePair{:});
end