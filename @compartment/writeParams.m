function writeParams(obj, filename, params, R, varargin)
% WRITEPARAMS is a private method for class multicompartment
%
%   writeParams store estimted parameters as a text file. 
%   Inpute:
%        obj: The multicompartment object
%   filename: Input filename to store results
%     params: A matrix of size #params X #voxels
%          R: An roi object. [optional]
%             Default: roi(#voxels, 1, 1)
%   
%   writeParams(..., <Parameter>, <Value>) allow controling the text format
%   by passing appropriate parameter-value pairs. 
%

% check inputs
[precision, mode] = parseInputVariables(varargin{:});
filename = validateFilename(filename);

% check if roi is provided
if exist('R', 'var')
    assert(isa(R, 'roi')&&length(R.index)==size(params,2), ...
        'MATLAB:multicompartment:writeParams',...
        'The thrid positional input argument must be an roi object.')
else
    R = roi(size(params,2), 1, 1);
end

%% write parameters
if strcmp(mode, 'new')
    if isfile(filename)
        warning('MATLAB:multicompartment:writeParams',...
                'Filename is overwritten.\n %s', filename);
    end
    fileId = fopen(filename,'w');
else
    fileId = fopen(filename,'a');
end

% write the opening header 
fprintf(fileId, '##parameters\n');

% write the column names
fprintf(fileId, '#$name 1 string\n');
colNames = ['nrow'; 'ncol'; 'ndep'; 'index'; obj.getParamsName];
fprintf(fileId, '%s ', colNames{:});
fprintf(fileId, '\n');

% write parameters for each voxel (one row per voxel)
nvoxels = size(params, 2);
fprintf(fileId, '#$value %d numeric\n', nvoxels);
for i=1:nvoxels
    fprintf(fileId, '%d %d %d %d', R.nrow, R.ncol, R.ndep, R.index(i));
    fprintf(fileId, sprintf(' %%1.%de', precision), params(:,i));
    fprintf(fileId, '\n');
end
fprintf(fileId, '##ENDparameters\n');
fclose(fileId);
end

function [precision, mode] = parseInputVariables(varargin)
p = inputParser;
p.CaseSensitive = false;
p.addParameter('precision'   , 6, @(v) assert(isnumeric(v)&&rem(v,1)==0));
p.addParameter('mode', 'new', @(v) assert(ischar(v)&&strcmp(v, 'new')||strcmp(v, 'append')));
p.parse(varargin{:});
precision = p.Results.precision;
mode = p.Results.mode;
end

function filename = validateFilename(filename)
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
end