function config = readConfig(filename, field)
    % readConfig is a method for class MULTICOMPARTMENT
    % 
    % readConfig(filename, field) read subfield for the given field and 
    % return a struct with subfield values in the format of cell.
    %             Input args:
    %               filename: Filename of type char to read the model
    %            Output args:
    %                  model: A struct
    %
    % Mohsen Farzi
    % Email: m.farzi@leeds.ac.uk
    
    % check arguments
    filename = validateFilename(filename);
    assert(isfile(filename), 'MATLAB:MULTICOMARTMENT:readConfig', ...
        'File does not exist:\n%s', filename);
    assert(ischar(field), 'MATLAB:MULTICOMARTMENT:readConfig', ...
        "Input field name must be of type 'char'.");
        
    % create output struct
    config = struct;
    thisHeader = strcat('##',field);
    headerExist = false;
    thisFooter = strcat('##END',field);
    footerExist = false;
    
    fileId = fopen(filename, 'r');
    thisLine = fgetl(fileId);
    while ischar(thisLine)
        if startsWith(thisLine, thisHeader) && ~headerExist
            headerExist = true;
        end
        
        if startsWith(thisLine, thisFooter) 
            footerExist = true;
            break;
        end
        
        if startsWith(thisLine, '#$') && headerExist
            [fieldname, fieldvalue, fileId] = readsubfield(fileId, thisLine);
            config.(fieldname) = fieldvalue;
        end
        
        thisLine = fgetl(fileId);
    end
    fclose(fileId);
    assert(headerExist, 'MATLAB:MULTICOMPARTMENT:readConfig',...
        'The header is missing for model.');
    assert(footerExist, 'MATLAB:MULTICOMPARTMENT:readConfig',...
        'The END is missing for model.');
end
    
function [fieldname, fieldvalue, fileId] = readsubfield(fileId, hdr)
    % read header
    c = textscan(hdr, '#$%s %d %s');
    fieldname = c{1}{1};
    assert(ischar(fieldname)&&~isempty(fieldname),...
        'MATLAB:MULTICOMPARTMENT:readsubfield',...
        "Corrupted header: field name must be of type 'char'.");
    N = c{2};
    assert(isnumeric(N) && isscalar(N) && N>0 && mod(N,1)==0,...
        'MATLAB:MULTICOMPARTMENT:readsubfield',...
        'Corrupted header: the number of lines must a positive integer.');
    type = c{3}{1};
    assert(ischar(type)&&(strcmp(type, 'string')||strcmp(type, 'numeric')),...
        'MATLAB:MULTICOMPARTMENT:readsubfield',...
        "Corrupted header: format must be either 'string' or 'numberic'.");
    if strcmp(type, 'string')
            format = '%s';
        else
            format = '%f';
    end
        
    % read field value
    fieldvalue = cell(N,1);
    for n=1:N
        thisLine = fgetl(fileId);
        c = textscan(thisLine, format);
        fieldvalue{n} = c{1};
    end
end