function add(obj, varargin)
    % add is a public method for class SCHEME
    %
    % add(obj, ghat, bval) add measuremetns for a scheme object of type
    % bvalue suitable for tensor fitting.
    %
    % add(obj, ghat, bval, dt, delta, te) add measurements for a scheme
    % object of type stejskal-tanner suitable for compartment modelling.
    %
    % add(...,<params>,<value>) allow setting nominal parameters.
    %
    %   input arguments:
    %              ghat: A double matrix [M x 3]. Columns represent the 
    %                    x,y, and z coordiantes and rows represent
    %                    diffusion measuremtns.
    %              bval: A non-negative double vector of size M. [s/m^2]
    %                dt: diffusion time [s]. Either a scalar or a postive
    %                    double vector of size M.
    %             delta: gradient duration [s]. Either a scaler or a 
    %                    postive double vector of size M.
    %                te: echo time [s]. Either a scaler or a positive 
    %                    double vector of size M.
    %
    %  input parameters:
    %       ghatNominal: A double matrix of size ghat.
    %       gmagNominal: A double vector of size bval. 
    

    % parse inputs
    p = inputParser;
    p.CaseSensitive = false;
    p.addRequired('ghat', @(v) isDirectionMatrix(v));
    p.addRequired('bval', @(v) isSemiPositiveVector(v));
    p.addOptional('dt'    , [], @(v) isPositiveVector(v));
    p.addOptional('delta' , [], @(v) isPositiveVector(v));
    p.addOptional('te'    , [], @(v) isPositiveVector(v));
    p.addParameter('ghatNominal', [], @(v) isDirectionMatrix(v));
    p.addParameter('gmagNominal', [], @(v) isSemiPositiveVector(v));
    p.addParameter('bvalNominal', [], @(v) isSemiPositiveVector(v));
    p.parse(varargin{:});
        
    ghat = scheme.normaliseDirections(p.Results.ghat);
    bval = p.Results.bval;
    dt   = p.Results.dt;
    delta= p.Results.delta;
    te   = p.Results.te;
    ghatNominal = p.Results.ghatNominal;
    gmagNominal = p.Results.gmagNominal;
    bvalNominal = p.Results.bvalNominal;
    
    nMeasurements = size(ghat, 1);    
    assert(size(bval,1)==nMeasurements,...
        'MATLAB:scheme:invalidInputArgument',...
        'Input ghat and bval must have equal number of rows.');
    
    optionalInputsGiven = [~isempty(dt), ~isempty(delta), ~isempty(te)];
    if strcmp(obj.type, 'bvalue')
        assert(all(~optionalInputsGiven), 'MATLAB:scheme:nargin',...
            'Too many input arguments for scheme type bvalue.');
        % update object 
        obj.ghat = [obj.ghat; ghat];
        obj.bval = [obj.bval; bval];
        
        % update nominal bvalue if provided
        if ~isempty(gmagNominal)
            error('MATLAB:scheme:invalidInputArg',...
                ['Set "bvalNominal" instead. "gmagNominal" is not',...
                ' supported for scheme type bvector']);
        elseif ~isempty(bvalNominal)
            assert(length(bvalNominal)==nMeasurements,...
            'MATLAB:scheme:invalidInputArg',...
            ['The length of bvalNominal must be equal to',...
             ' the number of rows in ghat.']);
            % update bvalDic
            obj.bvalDic = union(obj.bvalDic, bvalNominal, 'stable');
            [~, bvalCode] = ismember(bvalNominal,obj.bvalDic);
            assert(all(bvalCode>0),...
                'MATLAB:scheme:unknownError',...
                'Unknown error in clustering nominal bvalues.');
        else
            bvalCode = zeros(nMeasurements,1);
        end
        obj.bvalCode = [obj.bvalCode; bvalCode];
    else
        %type: stejskal-tanner
        assert(all(optionalInputsGiven), 'MATLAB:scheme:nargin',...
            'Too few input arguments for scheme type stejskal-tanner.');
        
        % check size for dt, delta, and te
        if isscalar(dt)
            dt = onese(nMeasurements, 1)*dt;
        else
            assert(length(dt)==nMeasurements,...
                'MATLAB:scheme:invalidInputArg',...
                ['Input dt must be either scalar or',...
                 ' its legnth equals the number of rows in ghat.']);
        end
        %
        if isscalar(delta)
            delta = onese(nMeasurements, 1)*delta;
        else
            assert(length(delta)==nMeasurements,...
                'MATLAB:scheme:invalidInputArg',...
                ['Input delta must be either scalar or',...
                 ' its legnth equals the number of rows in ghat.']);
        end
        %
        if isscalar(delta)
            delta = onese(nMeasurements, 1)*delta;
        else
            assert(length(delta)==nMeasurements,...
                'MATLAB:scheme:invalidInputArg',...
                ['Input te must be either scalar or',...
                 ' its legnth equals the number of rows in ghat.']);
        end

        % update scheme
        obj.ghat = [obj.ghat; ghat];
        
        obj.bval = [obj.bval; bval];
        
        obj.addDiffusionTimes(dt);
        [~,dtCode] = min(abs(dt-obj.dtList'),[],2);
        obj.dtCode = [obj.dtCode; dtCode];
        
        obj.addGradientDurations(delta);
        [~,deltaCode] = min(abs(delta-obj.deltaList'),[],2);
        obj.deltaCode = [obj.deltaCode; deltaCode];
        
        obj.addEchoTimes(te);
        [~,teCode] = min(abs(te-obj.teList'),[],2);
        obj.teCode = [obj.teCode; teCode];
        
        gmag = scheme.computeGmag(bval, dt, delta);
        obj.gmag = [obj.gmag; gmag];
        
        % update bvalCode and bvalDic 
        if ~isempty(gmagNominal)
            assert(isempty(bvalNominal),...
                'MATLAB:scheme:invalidInputArg',...
                'Either "bvalNominal" or "gmagNominal" must be provided.');
            assert(length(gmagNominal)==nMeasurements,...
                'MATLAB:scheme:invalidInputArg',...
                ['The length of gmagNominal must be equal to',...
                ' the number of rows in ghat.']);
            bvalNominal = scheme.computeBvalue(gmagNominal, dt, delta);
            % update bvalDic
            obj.bvalDic = union(obj.bvalDic, bvalNominal, 'stable');
            [~, bvalCode] = ismember(bvalNominal,obj.bvalDic);
            assert(all(bvalCode>0),...
                'MATLAB:scheme:unknownError',...
                'Unknown error in clustering nominal bvalues.');
        elseif ~isempty(bvalNominal)
            assert(length(bvalNominal)==nMeasurements,...
            'MATLAB:scheme:invalidInputArg',...
            ['The length of bvalNominal must be equal to',...
             ' the number of rows in ghat.']);
            % update bvalDic
            obj.bvalDic = union(obj.bvalDic, bvalNominal, 'stable');
            [~, bvalCode] = ismember(bvalNominal,obj.bvalDic);
            assert(all(bvalCode>0),...
                'MATLAB:scheme:unknownError',...
                'Unknown error in clustering nominal bvalues.');
        else
            bvalCode = zeros(nMeasurements,1);
        end
        obj.bvalCode = [obj.bvalCode; bvalCode];
    end
    
    % update nominal directions
    if ~isempty(ghatNominal)
        assert(size(ghatNominal,1)==nMeasurements,...
                'MATLAB:scheme:invalidInputArg',...
                'ghat and ghatNominal must have equal rows.');
                
        obj.ghatDic = union(obj.ghatDic, ghatNominal, 'stable', 'rows');
        [~, ghatCode] = ismember(ghatNominal, obj.ghatDic, 'rows');
        assert(all(ghatCode>0),...
             'MATLAB:scheme:unknownError',...
             'Unknown error in clustering nominal bvalues.');
    else
        ghatCode = zeros(nMeasurements, 1);
    end
    obj.ghatCode = [obj.ghatCode; ghatCode];
end





