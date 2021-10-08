classdef myo
    % MYO is a static class to read and write in '.myo' format
    
    methods (Static)
        config = read(filename, field);
        write(filename, config);
        filename = isValidFilename(filename);
    end
end