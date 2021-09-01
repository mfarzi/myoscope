function p = readParams(filename)
[filepath,name,ext] = fileparts(filename);
if isempty(filepath)
    filepath = pwd;
    filename = fullfile(filepath, strcat(name, ext));
end
if ~isfile(filename)
    error('MATLAB:compartment:readParams',...
          'Filename does not exist.\n%s', filename);
end
p = struct;
fileId = fopen(filename, 'r');
thisLine = fgetl(fileId);
while ischar(thisLine)
    if startsWith(thisLine, '#$')
        paramName = thisLine(3:end);
        paramVal = str2double(fgetl(fileId));
        p = setfield(p, paramName, paramVal);
    end

    thisLine = fgetl(fileId);
end
fclose(fileId);
end%of readParams