classdef myoscope 
    % MYOSCOPE
    %
    %   A static class to encapsulate methods to develope compartment
    %   models, fitting them to measured signals, read or write parameters,
    %   or visualise parameter maps.
    %
    %   methods (Static, public):
    %       getIndex          - get index 
    %       getVoxel          - return coordinates for voxels in the ROI
    %       update            - update ROI voxels
    %
    % Mohsen Farzi
    % Email: m.farzi@leeds.ac.uk
    
    methods (Static, Access='public') 
        model = str2model(modelName, params);
        model = load(fileName);
        save(model, fileName);
        property = searchConfig(fileName, varargin);
        %\\
    end%of methods
    %\\
end%of class MYOSCOPE