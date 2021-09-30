function filename = isValidFilename(filename)
    % check if input filename is a valid path to write the scheme object
    % into it
    
    % input type check
    assert(isa(filename, 'char'), 'MATLAB:myoscope:util',...
        'Input filename must be of type char.');
    
    % split filename 
    [filepath,name,ext] = fileparts(filename);
    % check folder exist
    if isempty(filepath)
        filepath = pwd;
    else
        assert(isfolder(filepath), 'MATLAB:myoscope:util',...
            'Input folder path does not exist.');
    end
    % check name is valid
    assert(isempty(regexp(name, '[^a-zA-Z_0-9\- ]', 'once')),...
        'MATLAB:myoscope:util', ['Input filename should only contain',...
        ' grouping of alphabetic, numeric, underscore, dash or space',...
        ' characters.']);
    % check extension
    if isempty(ext)
        ext = '.myo';
    elseif ~strcmp(ext, '.myo')
        ext = '.myo';
        warning('MATLAB:myoscope:util',...
            'The filename extension is changed to ".myo".');
    end
    
    % export filename
    filename = fullfile(filepath, strcat(name, ext));
end