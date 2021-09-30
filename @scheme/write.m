function write(obj, filename, varargin)
    % write(obj, filename, varargin) is a method for class SCHEME 
    
    M = obj.measurementsNum('all');
    assert(M>0, 'MATLAB:scheme:emptyObject',...
        'The scheme file must have at least one measurement.');
    
    % check if filename is valid
    filename = isValidFilename(filename);
    fileId = fopen(filename, 'w');
    
    % write scheme
    fprintf(fileId, '##scheme\n');
 
    fprintf(fileId, '#$type 1 string\n');
    fprintf(fileId, '%s\n', obj.type);   
    
    fprintf(fileId, '#$nMeasurements 1 numeric\n');
    fprintf(fileId, '%d\n', M);
    
    fprintf(fileId, '#$ghat 1 numeric\n');
    fmt = strcat(repmat('%1.6f ', 1, 3*M-1),' %1.6f\n');
    fprintf(fileId, fmt, obj.ghat');
    
    if strcmp(obj.type, 'bvalue')
        fprintf(fileId, '#$bval 1 numeric\n');
        fmt = strcat(repmat('%1.6e ', 1, M-1),' %1.6e\n'); 
        fprintf(fileId, fmt, obj.bval);
    else
        %stejskal-tanner
        fprintf(fileId, '#$gmag 1 numeric\n');
        fmt = strcat(repmat('%1.6e ', 1, M-1),' %1.6e\n'); 
        fprintf(fileId, fmt, obj.gmag);
        
        fprintf(fileId, '#$dt 2 numeric\n');
        N = length(obj.diffusionTimes);
        fmt = strcat(repmat('%1.6e ', 1, N-1), ' %1.6e\n');
        fprintf(fileId, fmt, obj.diffusionTimes);
        fmt = strcat(repmat('%d ', 1, M-1),' %d\n'); 
        fprintf(fileId, fmt, obj.dtCode);

        fprintf(fileId, '#$delta 2 numeric\n');
        N = length(obj.gradientDurations);
        fmt = strcat(repmat('%1.6e ', 1, N-1), ' %1.6e\n');
        fprintf(fileId, fmt, obj.gradientDurations);
        fmt = strcat(repmat('%d ', 1, M-1),' %d\n'); 
        fprintf(fileId, fmt, obj.deltaCode);


        fprintf(fileId, '#$te 2 numeric\n');
        N = length(obj.echoTimes);
        fmt = strcat(repmat('%1.6e ', 1, N-1), ' %1.6e\n');
        fprintf(fileId, fmt, obj.echoTimes);
        fmt = strcat(repmat('%d ', 1, M-1),' %d\n'); 
        fprintf(fileId, fmt, obj.teCode);
    end%if strcmp(obj.type, 'bvalue')
    
    % optional
    if ~isempty(obj.bvalDic)
        fprintf(fileId, '#$nominalBvals 2 numeric\n');
        nShells = size(obj.bvalDic, 1);
        fmt = strcat(repmat('%1.6e ', 1, nShells-1),' %1.6e\n'); 
        fprintf(fileId, fmt, obj.bvalDic);
        fmt = strcat(repmat('%d ', 1, M-1),' %d\n'); 
        fprintf(fileId, fmt, obj.bvalCode);
    end
    
    % optional
    if ~isempty(obj.ghatDic)
        fprintf(fileId, '#$nominalGhat 2 numeric\n');
        nD = size(obj.ghatDic, 1);
        fmt = strcat(repmat('%1.6f ', 1, 3*nD-1),' %1.6f\n'); 
        fprintf(fileId, fmt, obj.ghatDic');
        fmt = strcat(repmat('%d ', 1, M-1),' %d\n'); 
        fprintf(fileId, fmt, obj.ghatCode);
    end
    %
    fprintf(fileId, '##ENDscheme\n');
    fclose(fileId);
end