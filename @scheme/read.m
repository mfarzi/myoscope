function obj = read(filename)
    % read is a static method for class SCHEME
    % 
    % read(filename) return a SCHEME object.
    %
    % Mohsen Farzi
    % Email: m.farzi@leeds.ac.uk
    
    filename = isValidFilename(filename);
    assert(isfile(filename), 'MATLAB:scheme:read',...
        'Filename does not exist.');
    
    config = readConfig(filename, 'scheme');
    
    if strcmp(config.type{1}, 'bvalue')
        obj = scheme('bvalue');
        ghat = reshape(config.ghat{1}, 3, []);
        bval = config.bval{1};
        obj.add(ghat', bval);
    elseif strcmp(config.type{1}, 'stejskal-tanner')
        obj = scheme('stejskal-tanner');
        ghat = reshape(config.ghat{1}, 3, []);
        gmag = config.gmag{1};
        
        dtList = config.dt{1};
        dtCode = config.dt{2};
        dt = dtList(dtCode);
        
        deltaList = config.delta{1};
        deltaCode = config.delta{2};
        delta = deltaList(deltaCode);
        
        teList = config.te{1};
        teCode = config.te{2};
        te = teList(teCode);
            
        obj.add(ghat', gmag, dt, delta, te);
    else
        error('MATLAB:scheme:read',....
            'Unrecognised scheme type %s.', config.type{1});
    end
    
    % optional arguments
    if isfield(config, 'nominalBvals')
        obj.setNominalBvals(config.nominalBvals{1});
        assert(all(obj.bvalCode==config.nominalBvals{2}),...
            'MATLAB:scheme:inputFormat',...
            'Inconsistency in the input scheme file: nominalBvals.');
    end
    
    if isfield(config,'nominalGhat')
        ghatDic = reshape(config.nominalGhat{1}, 3, []);
        obj.setNominalDirections(ghatDic');
        assert(all(obj.ghatCode==config.nominalGhat{2}),...
            'MATLAB:scheme:inputFormat',...
            'Inconsistency in the input scheme file: nominalGhat.');
    end
end