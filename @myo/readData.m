function [sig, schemefile, roi] = readData(filename)
    filename = myo.isValidFilename(filename);

    if ~isfile(filename)
        error('MATLAB:myo:invalidFile',...
              'Filename does not exist.\n%s', filename);
    end
    
    config = myo.read(filename, 'data');
    sig = cell2mat(config.sig');
    nD = sig(1,1);
    roi = [sig(1:nD+1,1); sig(nD+2,:)'];
    sig = sig(nD+3:end,:);
    schemefile = scheme.read(filename);
end%of readParams