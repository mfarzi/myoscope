function newobj = subset(obj, idx)
    % subset return an independent copy of object with selected indices
    newobj = scheme(obj.type);
    ghat = obj.ghat(idx,:);
    bval = obj.bval(idx);
    dt = obj.dt(idx);
    delta = obj.delta(idx);
    te = obj.te(idx);
    
    ghatNominal = obj.ghatNominal(idx);
    bvalNominal = obj.bvalNominal(idx);
    
    newobj.add(ghat, bval, dt, delta, te, 'ghatNominal', ghatNominal,...
        'bvalNominal', bvalNominal);
end
    