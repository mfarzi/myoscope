function model = myoscopeLoadModel(fileName)
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

% read config file
config = searchConfig(fileName, 'name', 'constraints', 'params');

% construct the model
model = str2model(config.name);

defaultConstraints = model.getConstraints;

idx = ismember(config.constraints, defaultConstraints);
thisConstraints = config.constraints(~idx);
% enforce constraints
nConstraints = length(thisConstraints);
for i=1:nConstraints
    model.addConstraint(thisConstraints{i});
end

% set parameters values
model.set('params', config.params);
end % of loadCompartmentModel


    

