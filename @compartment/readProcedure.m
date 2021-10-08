function [model, sig, schemeFile] = readProcedure(filename)
    % readProcedure is a static method for the abstract class COMPARTMENT
    % 
    % readProcedure(filename) return the compartmen model, data, and
    % schemefile required to estimate parameters.
    %       Input args:
    %         filename: Filename of type char to read the procedure.
    %
    %      Output Args:
    %            model: A compartment object
    %              sig: A numeric matrix [nScheme x nVoxels]
    %       schemeFile: A scheme object with nScheme measurements
    %       
    %
    % see also: readProcedure, runProcedure
    
    filename = isValidFilename(filename);
    assert(isfile(filename), 'MATLAB:compartment:invalidInputArgument',...
        'Input filename does not exists.');
    
    model = compartment.readModel(filename);
    
    config = readConfig(filename, 'procedure');
    load(config.data{1}, 'sig');
    roi = config.data{2};
    
    schemeFile = scheme.read(config.scheme{1});
    idx = schemeFile.select(config.scheme{2});
    schemeFile.remove(not(idx));
    sig = sig(idx,:);
end

