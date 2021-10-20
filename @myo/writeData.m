function writeData(filename, sig, varargin)
    % WRITEPARAMS is a public method for class multicompartment
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
    
    

    % parse inputs
    p = inputParser;
    p.CaseSensitive = false;
    p.addOptional('schemefile', [], @(v) isa(v, 'scheme'));
    p.addOptional('roi', [], @(v) isvector(v) && all(mod(v,1)==0));
    p.addParameter('mode', 'a', @(v) ischar(v) && ismember(v,{'a', 'w'}));
    p.addParameter('precision', 10, @(v) isscalar(v)&&isnumeric(v)&&(mod(v,1)==0)&&v>0);
    p.parse(varargin{:});
    schemefile = p.Results.schemefile;
    roi = p.Results.roi;
    writeMode = p.Results.mode;
    precision = p.Results.precision;
    
    % check inputs
    if isempty(schemefile)
        schemefile = scheme.read(filename);
        printScheme = false;
    else
        printScheme = true;
    end
    
    if isempty(roi)
        nVoxels = size(sig, 2);
        roi = [1 nVoxels 1:nVoxels]';
    else
        nD = roi(1);
        assert(length(roi)==size(sig,2)+nD+1,...
        'MATLAB:compartment:invalidInputArg',...
        'roi has %d voxels but params has %d columns.',length(roi)-nD-1, size(sig,2));
    end
    
    assert(size(sig, 1) == schemefile.measurementsNum(),...
        'MATLAB:compartment:invalidInputArg',...
        'Input data matrix must have %d rows to match schemefile.',...
        schemefile.measurementsNum());
    
    % check filename is valid
    filename = myo.isValidFilename(filename);
    if strcmp(writeMode, 'w')
        % create new file and overwrite previous file
        fileId = fopen(filename, 'w');
    else
        % append data to the existing file
        fileId = fopen(filename, 'a');
    end

    %% write parameters
    % write the opening header 
    fprintf(fileId, '##data\n');
    
    fprintf(fileId, '#$sig %d numeric\n', size(sig,2));
    fmt = append('%.', num2str(precision), 'g');
    nD = roi(1);
    for i=1:size(sig,2)
        myo.print(fileId, [roi(1:nD+1); roi(i+nD+1); sig(:,i)], fmt);
    end
    fprintf(fileId, '##ENDdata\n\n');
    fclose(fileId);
    
    if printScheme
        schemefile.write(filename, 'mode', 'a');
    end
end