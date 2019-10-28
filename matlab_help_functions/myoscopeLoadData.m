function sig = myoscopeLoadData(path2data, path2roi)
% LOADCOMPARTMENTMODEL load compartment model 
%   model = LOADCOMPARTMENTMODEL(filename)
%   returns a multi-compartment object with given configurations.
%
% see also compartment

% Mohsen Farzi
% Email: m.farzi@leeds.ac.uk

% check arguments
if ~isfile(path2data)
    error('MATLAB:loadCompartmentModel:fileIsNotFound', ...
          'File does not exist:\n%s', path2data);
end

if ~isfile(path2roi) && ~isempty(path2roi)
    error('MATLAB:loadCompartmentModel:fileIsNotFound', ...
          'File does not exist:\n%s', path2roi);
end

load(path2data, 'data');
if isempty(path2roi)
    sig = data;
else
    load(path2roi, 'roi');
    idx = roi(roi>0);
    sig = data(:, idx);
end
end % of loadCompartmentModel


    

