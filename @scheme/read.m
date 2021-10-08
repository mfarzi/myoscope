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
    
    %% read config
    type = config.type{1};
    assert(ismember(type, {'bvalue', 'stejskal-tanner'}),...
        'MATLAB:scheme:unknowFileFormat',...
        'Unrecognised scheme type %s.', type);
    
    nMeasurements = config.nMeasurements{1};
    
    ghat = reshape(config.ghat{1}, 3, []);
    bval = config.bval{1};
    
    if isfield(config, 'dt')
        dtList = config.dt{1};
        dtCode = config.dt{2};
        dt = dtList(dtCode);
    end
    
    if isfield(config, 'delta')
        deltaList = config.delta{1};
        deltaCode = config.delta{2};
        delta = deltaList(deltaCode);
    end
    
    if isfield(config, 'te')
        teList = config.te{1};
        teCode = config.te{2};
        te = teList(teCode);
    end    
    %% create scheme object
    % bvalue
    if strcmp(type, 'bvalue')
        obj = scheme('bvalue');
        obj.add(ghat', bval);
    else
        obj = scheme('stejskal-tanner');
        obj.add(ghat', bval, dt, delta, te);
    end
    
    if isfield(config, 'nominalBvals')
        obj.bvalDic = config.nominalBvals{1};
        obj.bvalCode = config.nominalBvals{2};
    else
        obj.bvalCode = zeros(nMeasurements, 1);
    end
    
    if isfield(config,'nominalGhat')
        ghatDic = reshape(config.nominalGhat{1}, 3, []);
        obj.ghatDic = ghatDic';
        obj.ghatCode = config.nominalGhat{2};
    else
        obj.ghatCode = zeros('nMeasurements', 1);
    end 
end