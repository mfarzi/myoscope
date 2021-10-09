function writeProcedure(obj, filename, path2data, path2scheme, roi, schemeFilter)
    % writeProcedure is a method for the abstract class COMPARTMENT
    % 
    % writeProcedure(obj,filename,path2data,path2scheme,[roi],[schemeFilter]) 
    % save the model configuration, data path, diffusion scheme path, and 
    % optional roi and scheme filter. Method 'runProcedure' may be employed to
    % fit model parameters to DW-MR data and save the parameters.
    %       Input args:
    %              obj: A compartment object
    %         filename: Filename of type char to write the procedure.
    %        path2data: Path to data in MATLAB format
    %      path2scheme: Path to scheme file
    %            [roi]: Optional roi object to define a region of interest
    %   [schemeFilter]: Optional filter of type char to select a subset of
    %                   desired diffusion measurements           
    %
    % see also: readProcedure, runProcedure

    % parse input arguments
    filename = myo.isValidFilename(filename);
    assert(~isfile(filename), 'MATLAB:compartment:invalidInputArgument',...
        'Input filename already exists.');

    %% initiate a fileId
    fileId = fopen(filename,'w');

    % write headers
    fprintf(fileId, '##filepath\n');

    % write data path
    fprintf(fileId, '#$data 1 string\n');
    fprintf(fileId, '%s\n', path2data);

    % write roi
    fprintf(fileId, '#$roi 1 numeric\n');
    myo.print(fileId, roi, '%d');
    
    % write scheme file path
    fprintf(fileId, '#$scheme 1 string\n');
    fprintf(fileId, '%s\n', path2scheme);
    
    % write schemfilter
    fprintf(fileId, '#$filter 1 numeric\n');
    myo.print(fileId, schemeFilter, '%d');

    fprintf(fileId, '##ENDfilepath\n\n');
    fclose(fileId);
    
    % add model config to the filename
    obj.writeModel(filename, 'a');
end

