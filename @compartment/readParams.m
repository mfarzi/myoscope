function [params, rmse, roi] = readParams(filename)
    filename = myo.isValidFilename(filename);

    if ~isfile(filename)
        error('MATLAB:compartment:readParams',...
              'Filename does not exist.\n%s', filename);
    end
    
    config = myo.read(filename, 'parameters');
    params = cell2mat(config.params');
    % read roi
    nD = params(1,1);
    roi = [params(1:nD+1,1); params(nD+2,:)'];
    params = params(nD+3:end,:);
    
    if isfield(config, 'rmse')
        rmse = cell2mat(config.rmse');
        rmse = rmse(nD+3,:)';
    else
        rmse = [];
    end
end%of readParams