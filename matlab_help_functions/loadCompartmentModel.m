function model = loadCompartmentModel(fileName)
% LOADCOMPARTMENTMODEL load compartment model 
%   model = LOADCOMPARTMENTMODEL(filename)
%   returns a multi-compartment object with given configurations.
%
% see also compartment

% Mohsen Farzi
% Email: m.farzi@leeds.ac.uk

% check arguments
if ~isfile(fileName)
    error('MATLAB:loadCompartmentModel:fileIsNotFound', ...
          'File does not exist:\n%s', fileName);
end

[config, flag] = searchConfig(fileName, 'name');
model = str2model(config.name{1});
end % of loadCompartmentModel


    

