function writeParams(obj, filename, params, varargin)
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
    p.addOptional('rmse', [], @(v) isPositiveVector(v));
    p.addOptional('roi', [], @(v) isvector(v) && all(mod(v,1)==0));
    p.addParameter('mode', 'a', @(v) ischar(v) && ismember(v,{'a', 'w'}));
    p.addParameter('precision', 10, @(v) isscalar(v)&&isnumeric(v)&&(mod(v,1)==0)&&v>0);
    p.parse(varargin{:});
    rmse = p.Results.rmse;
    roi = p.Results.roi;
    writeMode = p.Results.mode;
    precision = p.Results.precision;
    
    % check inputs
    assert(isempty(rmse)||length(rmse)==size(params,2),...
        'MATLAB:compartment:invalidInputArg',...
        'rmse and params must have eqaul number of columns.');
    
    if isempty(roi)
        nVoxels = size(params, 2);
        roi = [1 nVoxels 1:nVoxels]';
    else
        nD = roi(1);
        assert(length(roi)==size(params,2)+nD+1,...
        'MATLAB:compartment:invalidInputArg',...
        'roi and params must have eqaul number of columns.');
    end
    
    assert(size(params, 1) == obj.getParamsNum(),...
        'MATLAB:compartment:invalidInputArg',...
        ['Matrix params must have %d rows including ',...
         'roi information.'], obj.getParamsNum());
    
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
    fprintf(fileId, '##parameters\n');
    
    fprintf(fileId, '#$params %d numeric\n', size(params,2));
    fmt = append('%.', num2str(precision), 'g');
    nD = roi(1);
    for i=1:size(params,2)
        myo.print(fileId, [roi(1:nD+1); roi(i+nD+1); params(:,i)], fmt);
    end
    
    if ~isempty(rmse)
        fprintf(fileId, '#$rmse %d numeric\n', size(params, 2));
        for i=1:size(params,2)
            myo.print(fileId, [roi(1:nD+1); roi(i+nD+1); rmse(i)], fmt);
        end
    end
    fprintf(fileId, '##ENDparameters\n\n');
    fclose(fileId);
end