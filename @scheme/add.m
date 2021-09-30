function add(obj, varargin)
    % add is a public method for class SCHEME
    %
    % add(obj, ghat, bval) return a scheme object of type bvalue suitable
    % for tensor fitting.
    %
    % add(obj, ghat, gmag, dt, delta, [te]) return a scheme object of type
    % stejskal-tanner.
    %
    %   input arguments:
    %              ghat: A double matrix. Columns represent the x,y, and z
    %                    coordiantes and rows represent diffusion
    %                    measuremtns.
    %              bval: A positive double vector of size M. [s/m^2]
    %              gmag: A positive double vector of size M. [T/m]
    %                dt: diffusion time [s]. A postive double vector of 
    %                    size M.
    %             delta: gradient duration [s]. A postive double vector of
    %                    size M.
    %                te: echo time [s]. A positive double vector of size M.
    %
    %% return type: bvalue
    if strcmp(obj.type, 'bvalue')
        p = inputParser;
        p.CaseSensitive = false;
        p.addRequired('ghat', @(v) isDirectionMatrix(v));
        p.addRequired('bval', @(v) isSemiPositiveVector(v));
        p.parse(varargin{:});
        
        ghat = scheme.normaliseDirections(p.Results.ghat);
        bval = p.Results.bval;
        
        assert(size(ghat,1)==size(bval,1),...
            'MATLAB:scheme:addMeasurements',...
            'Input ghat and bval must have equal number of rows.');
        
        % update object 
        obj.ghat = [obj.ghat; ghat];
        obj.bval = [obj.bval; bval];
    end
    
    %% return type: stejskal-tanner
    if strcmp(obj.type, 'stejskal-tanner')
        p = inputParser;
        p.CaseSensitive = false;
        p.addRequired('ghat'  , @(v) isDirectionMatrix(v));
        p.addRequired('gmag'  , @(v) isSemiPositiveVector(v));
        p.addRequired('dt'    , @(v) isPositiveVector(v));
        p.addRequired('delta' , @(v) isPositiveVector(v));
        p.addRequired('te'    , @(v) isPositiveVector(v));
        p.parse(varargin{:});
        
        ghat = scheme.normaliseDirections(p.Results.ghat);
        gmag = p.Results.gmag;
        dt   = p.Results.dt;
        delta= p.Results.delta;
        te   = p.Results.te;
        bval = scheme.computeBvalue(gmag, dt, delta);
        
        nrows = size(ghat, 1);
        assert(size(gmag,1)==nrows && size(delta,1)==nrows && ...
            size(dt,1)==nrows && size(te,1)==nrows,...
            'MATLAB:scheme:addMeasurements',...
            'Input arguments must have equal number of rows.');
        
        % update scheme
        obj.ghat = [obj.ghat; ghat];
        obj.gmag = [obj.gmag; gmag];
        obj.bval = [obj.bval; bval];
        
        obj.addDiffusionTimes(dt);
        obj.addGradientDurations(delta);
        obj.addEchoTimes(te);
        
        [~,dtCode] = min(abs(dt-obj.diffusionTimes'),[],2);
        obj.dtCode = [obj.dtCode; dtCode];
        
        [~,deltaCode] = min(abs(delta-obj.gradientDurations'),[],2);
        obj.deltaCode = [obj.deltaCode; deltaCode];
        
        [~,teCode] = min(abs(te-obj.echoTimes'),[],2);
        obj.teCode = [obj.teCode; teCode];
    end
end





