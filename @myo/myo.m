classdef myo
    % MYO is a static class to read and write in '.myo' format
    
    methods (Static)
        config = read(filename, field);
        write(fielId, fieldname, v);
        filename = isValidFilename(filename);
        print(fileId, v, fmt0);
        writeData(filename, sig, varargin);
        [sig, schemefile, roi] = readData(filename);
    end
end