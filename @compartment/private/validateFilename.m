function filename = validateFilename(filename)
% validateFilename is a private method for class MULTICOMPARTMENT
%
% validateFielanme(filename) asserts the validity of the input attributes.

assert(isa(filename, 'char'),...
    'MATLAB:MULTICOMPARTMENT:validateFilename',...
    'Input filename must be of type char.');

[filepath,name,ext] = fileparts(filename);
if isempty(filepath)
    filepath = pwd;
end
if ~strcmp(ext, '.myo')
    error('MATLAB:MULTICOMPARTMENT:validateFilename',...
         "'The filename extension must be '.myo'.");    
end
filename = fullfile(filepath, strcat(name,ext));