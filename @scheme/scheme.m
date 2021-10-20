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
    %       gmag        - gradient magnitude [T/m]
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
        
        dtList = [];         % Difussion times [s]
        dtCode = [];         % 
        
        deltaList = [];      % gradient durations [s]
        deltaCode = [];
        
        teList = [];         % echo times [s]
        teCode = [];
        
        ghatDic = zeros(0,3);
        ghatCode = [];
        
        bvalDic = [];
        bvalCode = [];
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
        
        function v = dt(obj, idx)
            if nargin==1
                idx = true(obj.measurementsNum, 1);
            end
            
            if strcmp(obj.type, 'bvalue')
                v = [];
            else
                v = obj.dtList(obj.dtCode(idx));
            end
        end
        
        function v = delta(obj, idx)
            if nargin==1
                idx = true(obj.measurementsNum, 1);
            end
            
            if strcmp(obj.type, 'bvalue')
                v = [];
            else
                v = obj.deltaList(obj.deltaCode(idx));
            end
        end
        
        function v = te(obj, idx)
            if nargin==1
                idx = true(obj.measurementsNum, 1);
            end
            
            if strcmp(obj.type, 'bvalue')
                v = [];
            else
                v = obj.teList(obj.teCode(idx));
            end
        end
        
        function g = ghatNominal(obj, idx)
            if nargin==1
                idx = true(obj.measurementsNum, 1);
            end
            
            if any(obj.ghatCode==0)
                g = zeros(0,3);
            else
                g = obj.ghatDic(obj.ghatCode(idx),:);
            end
        end
        
        function b = bvalNominal(obj, idx)
            if nargin==1
                idx = true(obj.measurementsNum, 1);
            end
            
            if any(obj.bvalCode==0)
                b = [];
            else
                b = obj.bvalDic(obj.bvalCode(idx));
            end 
        end
        
        n = measurementsNum(obj,str);
        write(obj, filename, varargin);
        tbl = exportAsTable(obj);
        add(obj, varargin);
        setNominalBvals(obj, b);
        setNominalDirections(obj, ghat);
        idx = select(obj, str);
        out = get(obj, str);
        remove(obj, idx);
        newobj = subset(obj, idx);
        obj = vertcat(obj, schemefile);
    end
    
    methods (Static)
        [obj, gVec] = read(filename);
        [ghatCode, ghatDic] = clusterDirections(ghat, varargin);
        g = normaliseDirections(ghat);
        b = computeBvalue(gmag, dt, delta);
        g = computeGmag(bval, dt, delta);
    end
end%of class scheme
    