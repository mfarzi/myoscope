function write(obj, filename, varargin)
    % write(obj, filename, varargin) is a method for class SCHEME 
    
    % parse inputs
    p = inputParser;
    p.CaseSensitive = false;
    p.addOptional('writePermission', 'w', @(v) ischar(v) && ismember(v,{'a', 'w'}));
    p.addParameter('precision', 10, @(v) isscalar(v)&&isnumeric(v)&&(mod(v,1)==0)&&v>0);
    p.parse(varargin{:});
    
    writePermission = p.Results.writePermission;
    precision = p.Results.precision;
    
    M = obj.measurementsNum('all');
    assert(M>0, 'MATLAB:scheme:emptyObject',...
        'The scheme file must have at least one measurement.');
    
    % check if filename is valid
    filename = isValidFilename(filename);
    if strcmp(writePermission, 'w')
        % create new file and overwrite previous file
        fileId = fopen(filename, 'w');
    else
        % append data to the existing file
        fileId = fopen(filename, 'a');
    end
    
    % write scheme
    fprintf(fileId, '##scheme\n');
 
    fprintf(fileId, '#$type 1 string\n');
    fprintf(fileId, '%s\n', obj.type);   
    
    fprintf(fileId, '#$nMeasurements 1 numeric\n');
    fprintf(fileId, '%d\n', M);
    
    fprintf(fileId, '#$ghat 1 numeric\n');
    str = myprint(obj.ghat', precision);
    fprintf(fileId, '%s', str);
    
    fprintf(fileId, '#$bval 1 numeric\n');
    str = myprint(obj.bval, precision);
    fprintf(fileId, '%s', str);
    
    if strcmp(obj.type, 'stejskal-tanner')
        % write diffusion times, gradient durations, and echo times
        fprintf(fileId, '#$dt 2 numeric\n');
        str = myprint(obj.dtList, precision);
        fprintf(fileId, '%s', str);
        str = myprint(obj.dtCode, precision);
        fprintf(fileId, '%s', str);

        fprintf(fileId, '#$delta 2 numeric\n');
        str = myprint(obj.deltaList, precision);
        fprintf(fileId, '%s', str);
        str = myprint(obj.deltaCode, precision);
        fprintf(fileId, '%s', str);

        fprintf(fileId, '#$te 2 numeric\n');
        str = myprint(obj.teList, precision);
        fprintf(fileId, '%s', str);
        str = myprint(obj.teCode, precision);
        fprintf(fileId, '%s', str);
    end
    
    % optional inputs
    if ~isempty(obj.bvalDic)
        fprintf(fileId, '#$nominalBvals 2 numeric\n');
        str = myprint(obj.bvalDic, precision);
        fprintf(fileId, '%s', str);
        str = myprint(obj.bvalCode, precision);
        fprintf(fileId, '%s', str);
    end
    
    % optional inputs
    if ~isempty(obj.ghatDic)
        fprintf(fileId, '#$nominalGhat 2 numeric\n');
        str = myprint(obj.ghatDic', precision);
        fprintf(fileId, '%s', str);
        str = myprint(obj.ghatCode, precision);
        fprintf(fileId, '%s', str);
    end
    %
    fprintf(fileId, '##ENDscheme\n');
    fclose(fileId);
end