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
    p.addOptional('writePermission', 'a', @(v) ischar(v) && ismember(v,{'a', 'w'}));
    p.addParameter('precision', 10, @(v) isscalar(v)&&isnumeric(v)&&(mod(v,1)==0)&&v>0);
    p.parse(varargin{:});
    rmse = p.Results.rmse;
    writePermission = p.Results.writePermission;
    precision = p.Results.precision;
    
    % check inputs
    assert(isempty(rmse)||length(rmse)==size(params,2),...
        'MATLAB:compartment:invalidInputArg',...
        'rmse and params must have eqaul number of columns.');
    nD = params(1,1);
    assert(mod(nD,1)==0 && all(params(1,:)==nD) && nD>0,...
        'MATLAB:compartment:invalidInputArg',...
        ['The first row of matrix params must define ',...
         'the ROI dimension. For example, for 3D images, it must be 3.']);
    assert(size(params, 1) == obj.getParamsNum()+nD+2,...
        'MATLAB:compartment:invalidInputArg',...
        ['Matrix params must have %d rows including ',...
         'roi information.'], obj.getParamsNum()+nD+2);
    
    % check filename is valid
    filename = isValidFilename(filename);
    if strcmp(writePermission, 'w')
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
    for i=1:size(params,2)
        str = myprint(params(:,i), precision);
        fprintf(fileId, '%s', str);
    end
    
    if ~isempty(rmse)
        fprintf(fileId, '#$rmse %d numeric\n', size(params, 2));
        for i=1:size(params,2)
            str = myprint([params(1:5,i);rmse(i)], precision);
            fprintf(fileId, '%s', str);
        end
    end
    fprintf(fileId, '##ENDparameters\n\n');
    fclose(fileId);
end