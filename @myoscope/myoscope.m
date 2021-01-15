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
    
    properties (Access='private')
        model = [];         % mulitcompartment object
        cwd = pwd;
        params = [];
        rmse = [];
        aic = [];
        voxels = [];
    end
    
    methods
        function obj = myoscope(varargin)
            %MYOSCOPE Construct Function
            if nargin==0
                %Do Nothing
            else
                p = inputParser();
                p.CaseSensitive = false;
                p.addOptional('name', [], @(v) assert(ischar(v)));
                p.addParameter('cwd', pwd);
                
                p.parse(varargin{:});
                
                modelName = p.Results.name;
                obj.model = myoscope.str2model(modelName);
                
                cwd = p.Results.cwd;
                if isfolder(cwd)
                    obj.cwd = cwd;
                else
                    obj.cwd = pwd;
                    warning('MATLAB:MYOSCOPE:Construction',...
                            strcat('Input directory is not valid. The',...
                            ' current working directory is set as\n',...
                            '%s'),pwd);
                end
             end
        end
        writeModel(obj, filename, varargin);
        writeParams(obj, filename);
        obj = fit(obj, scheme, data, varargin);
    end
    
    %methods (Static, Access='public') 
    %    obj = readModel(fileName);
        %\\
    %end%of methods
    
    methods (Static, Access='private')
        model = str2model(modelName, params);
        property = searchConfig(fileName, varargin);
    end
    %\\
end%of class MYOSCOPE