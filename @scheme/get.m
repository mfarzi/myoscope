function out = get(obj, str)
    % return class properties
    assert(ischar(str)&&~isempty(str),...
        'MATLAB:scheme:invalidInputArgument',...
        'Input must be of type char.');
    
    switch str
        case 'ghatCode'
            out = obj.ghatCode;
        case 'bval'
            out = obj.bval;
        case 'dt'
            out = obj.dt;
        case 'te'
            out = obj.te;
        case 'ghat'
            out = obj.ghat;
        case 'delta'
            out = obj.delta;
        case 'ghatNominal'
            out = obj.ghatNominal;
        case 'ghatDic'
            out = obj.ghatDic;
        case 'bvalNominal'
            out = obj.bvalNominal;
        otherwise
            error('MATLAB:scheme:invalidInputArgument',...
                'Input %s is not valid.', str);
    end
end