classdef scheme < handle
    % SCHEME
    %
    %   A SCHEME object facilitate the management of different diffusion
    %   scheme, i.e. stejskal-tanner or b-vecotr. 
    %
    %   properties:
    %       type        - 'stejskal-tanner' or 'bvalue'
    %       ghat        - gradient directions [Mx3] numeric matrix
    %       bval        - diffusion coefficient; b-value [s/m^2]
    %       gmag        - gradient magnitude
    %
    %   methods (public):
    %       dt          - return diffusion times [s]  
    %       delta       - return gradient durations [s]
    %       te          - return echo times.
    %       bvalNominal - return nominal b-values. Due to gradient crushers
    %                     the experimental and expected nominal b-values
    %                     are not the same.
    %       ghatNominal - return nominal gradient directions.
    %       add         - Add new measurements to scheme object
    %       remove      - Remove selected measurements
    %       select      - Return indices to the selected measurements
    %       setNominalBvals - 
    %       write       - write scheme file in a text file.
    %       read        - read scheme file from an input text file.
    %       summary     - return a summary of the scheme object.
    %
    % Mohsen Farzi
    % Email: m.farzi@leeds.ac.uk
    
    properties (SetAccess = 'private')
        type = [];           % 'stejskal-tanner' or 'bvalue'
        ghat = [];           % Actual coordinates for the gradient
                             % directions; Double matrix of size Mx3 
        gmag = [];           % gradient magnitudes [T/m]
        bval = [];           % experimental diffusion coefficient [s/m^2]
        
        ghatDic = [];
        ghatCode = [];
        
        bvalDic = [];
        bvalCode = [];
    end
    
    properties (Access='private')
        diffusionTimes = [];
        gradientDurations = [];
        echoTimes = [];
        
        dtCode = [];
        deltaCode = [];
        teCode = [];
    end
        
    
    methods
        function obj = scheme(str)
            %SCHEME Construct Function.
            %
            %   scheme(bw) construct an empty object
            if nargin==0 || isempty(str)
                obj.type = 'stejskal-tanner';
            elseif ischar(str)&& ismember(str, {'stejskal-tanner', 'bvalue'})
                obj.type = str;
            else
                error('MATLAB:scheme:constructor',...
                    "Input type must be 'stejskal-tanner' or 'b-value'.");
            end
        end
        
        function v = dt(obj)
            v = obj.diffusionTimes(obj.dtCode);
        end
        
        function v = delta(obj)
            v = obj.gradientDurations(obj.deltaCode);
        end
        
        function v = te(obj)
            v = obj.echoTimes(obj.teCode);
        end
        
        function g = ghatNominal(obj)
            assert(~isempty(obj.ghatCode),...
                'MATLAB:scheme:unsetProperty',...
                ['Nominal gradient directions are not set yet. Use',...
                ' "setNominalDirections" method.']);
            g = obj.ghatDic(obj.ghatCode,:);
        end
        
        function b = bvalNominal(obj)
            assert(~isempty(obj.bvalCode),...
                'MATLAB:scheme:unsetProperty',...
                ['Nominal B-values are not set yet. Use',...
                ' "setNominalBvals" method.']);
            b = obj.bvalDic(obj.bvalCode);
        end
        
        n = measurementsNum(obj,str);
        write(obj, filename);
        tbl = exportAsTable(obj);
        add(obj, varargin);
        setNominalBvals(obj, b);
        setNominalDirections(obj, ghat);
        idx = select(obj, str);
        out = get(obj, str);
        remove(obj, idx);
    end
    
    methods (Static)
        [obj, gVec] = read(filename);
        [ghatCode, ghatDic] = clusterDirections(ghat, varargin);
        g = normaliseDirections(ghat);
        b = computeBvalue(gmag, dt, delta)
    end
end%of class scheme
    