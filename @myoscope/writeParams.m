function writeParams(obj, filename, varargin)
% WRITEPARAMS is a static method for class myoscope
%
% myoscope.writeParams(params, [roi], [model], [scheme], [fitter], [data]) 
% store parameters in a text file with extension 'myo'. Optional positional 
% arguments are allowed.
%           params: A matrix of size #params X #voxels         
%
%% parse input parameters
p = inputParser;
p.CaseSensitive = false;
p.addParameter('precision'   , 6, @(v) assert(isnum(v)&&rem(v,1)==0));
p.addOptional('mode', 'new', @(v) assert(strcmp(v, 'new')||strcmp(v, 'append')));
p.parse(varargin{:});
precision = p.Results.precision;
mode = p.Results.mode;

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
%% write parameters
if strcmp(mode, 'new')
    if isfile(filename)
        warning('MATLAB:MYOSCOPE:writeParams',...
                'Filename is overwritten.\n %s', filename);
    end
    fileId = fopen(filename,'w');
else
    fileId = fopen(filename,'a');
end

% write parameters
fprintf(fileId, '##parameters\n');
paramsList = obj.model.getParamsList;
paramsNum = length(paramsList);
for i=1:paramsNum
    fprintf(fileId, '#$%s\n', paramsList{i});
    fprintf(fileId, '%1.6e ', obj.params(:,i));
    fprintf(fileId, '\n');
end
fprintf(fileId, '##end\n');
fclose(fileId);
% write ROI
obj.voxels.write(filename,'append');
% write model
obj.writeModel(filename, 'append');
end