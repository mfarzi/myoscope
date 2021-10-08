function obj = vertcat(obj, schemefile)
    % a public method for class scheme
    assert(isa(schemefile, 'scheme'),...
            'MATLAB:scheme:invalidInputArgument',...
            'Cannot concatanate scheme class with %s.',...
            class(schemefile));
        
    assert(strcmp(obj.type, schemefile.type),...
            'MATLAB:scheme:invalidInputArgument',...
            'Input scheme files must have the same type.');
    
    % read parameters
    ghat = schemefile.ghat;
    bval = schemefile.bval;
    dt = schemefile.dt;
    delta = schemefile.delta;
    te = schemefile.te;
    
    ghatNominal = schemefile.ghatNominal;
    bvalNominal = schemefile.bvalNominal;
    
    obj.add(ghat, bval, dt, delta, te,...
            'ghatNominal', ghatNominal, ...
            'bvalNominal', bvalNominal);
end