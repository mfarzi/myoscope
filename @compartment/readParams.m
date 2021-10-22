function [params, rmse] = readParams(filename)
    filename = myo.isValidFilename(filename);

    if ~isfile(filename)
        error('MATLAB:compartment:readParams',...
              'Filename does not exist.\n%s', filename);
    end
    
    config = myo.read(filename, 'parameters');
    params = cell2mat(config.params');
    if isfield(config, 'rmse')
        rmse = cell2mat(config.rmse');
    else
        rmse = [];
    end
end%of readParams