function [params, rmse] = readParams(filename)
    filename = myo.isValidFilename(filename);

    if ~isfile(filename)
        error('MATLAB:compartment:readParams',...
              'Filename does not exist.\n%s', filename);
    end
    
    config = myo.read(filename, 'parameters');
    params = cell2mat(config.params');
    rmse = cell2mat(config.rmse');
end%of readParams