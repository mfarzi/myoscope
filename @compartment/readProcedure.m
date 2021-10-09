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
    
    filename = myo.isValidFilename(filename);
    assert(isfile(filename), 'MATLAB:compartment:invalidInputArgument',...
        'Input filename does not exists.');
    
    model = compartment.readModel(filename);
    
    config = readConfig(filename, 'filepath');
    load(config.data{1}{1}, 'sig');
    roi = config.roi{1};
    nD = roi(1);
    for i=1:nD+1
        assert(all(sig(i,:)==roi(i)),...
            'MATLAB:compartment:inconsistentFormat',...
            'Input ROI and singal formats are not consistent.');
    end
    iVoxels = ismember(sig(nD+2,:), roi(nD+2:end));
    sig = sig(:,iVoxels);
    
    thisScheme = scheme.read(config.scheme{1}{1});
    idx = config.filter{1}>0;
    schemeFile = thisScheme.subset(idx);
    
    idx = [true(nD+2,1); idx];
    sig = sig(idx,:);
end

