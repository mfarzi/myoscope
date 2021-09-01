classdef linker < matlab.mixin.Copyable
    % LINKER 
    %   
    %   A LINKER object allows mapping an unconstrained optimisation 
    %   variable x to a constrained model parameter p to be used in 
    %   conjuction with COMPARTMENT classes. Linker objects can be 
    %   concatenated in column arrays using vertcat to assist mapping 
    %   vector inputs X. 
    %   Note in this documantation, captial X and P represent column vector
    %   variables of size N whereas x and p represent scalar optimisation 
    %   variable and constrained model parameters, respectively.
    % 
    %
    %   properties (public):
    %       name           - semantic name given to variable x
    %       bounded        - Boolean
    %       type           - 'independent'
    %                           p[n] depends only on x[n]
    %                        'dependant'
    %                           p[n] depends on x[n] but its upperbound or
    %                           lowerbound depends on other parameters 
    %                        'dummy'
    %                           p[n] depends only on other parameters 
    %                           rather than x[n]
    %                        'constant'
    %                           p[n] is fixed at a given value.
    %       maximum        - max value for parameter p
    %       minimum        - min value for parameter p
    %       constraint     - constraint enforced on variable x
    %       dependence     - indices of variables on which x dependends
    %       compOrder      - computation order
    %
    %   methods (public):
    %       map            - forward transform p = F(x)  
    %       invmap         - backward transform x = F^(-1)(p)
    %       jacobian       - jacobian of P wrt x
    %
    % Mohsen Farzi
    % Email: m.farzi@leeds.ac.uk
    
    properties (Access = 'public')
        name = [];          % semantic variable name
        bounded = false;    % true or false  
        type ='independent';% 'independent': parameter p[i] depends only on
                            %                x[i]
                            %   'dependant': parameter p[i] depends on x[i]
                            %                and other parameters p[j],...
                            %       'dummy': parameter p[i] depends only on
                            %                other parameters p[j],...
                            %    'constant': parameter p[i] is fixed at a
                            %                given value.
        maximum = [];       % max value for 'bounded' variables
        minimum = [];       % min value for 'bounded' variables
        constraint = '';    % constraint enforced on variable x
        dependence = [];    % indices of variables on which x[n] dependends
        compOrder = 1;      % computation order of linking parameters;
                            % integer number from 1 to N (object length)  
    end
    
    properties (Access='public')   
        % function handles
        fval = [];          % @(x, P, u, l) map x to p 
        fvalinv = [];       % @(p, u, l) map p to x
        fgrad = [];         % @(x, u, l) dp/dx     
        
        ubound = [];        % @(P, u, l) return the upperbound
        uboundgrad = [];    % @(jac, P, N) return d(ubound)/dX
        lbound = [];        % @(P, u, l) return the lowerbound
        lboundgrad = [];    % @(jac, P, N) return d(ubound)/dX
        
        index = 1;          % index of this object in a column array
        varNum = 1;         % Length of the parrent object                     
    end
    
    methods
        function obj = linker(varargin)
            %LINKER Construct Function.
            %
            %   linker(name, 'unbounded') creates a linear link with given
            %   name and p=x identity mapping function.
            %
            %   linker(name, 'bounded', lowerbound, upperbound) creates a 
            %   link with given name and p = (u-l)*cos^2(x)+l mapping 
            %   function.
            %
            %   linker(..., <Parameter>, <Value>, ...) allow
            %   passing model parameters as parameter-value pairs.
            %        
            
            % parse input arguments
            args = parseInputArgs(varargin{:});
            
            obj.name = args.name;
            obj.bounded = strcmp(args.type, 'bounded');
            obj.type = 'independent';
            obj.maximum = args.maximum;
            obj.minimum = args.minimum;
            
            if obj.bounded
                obj.fval = @(x,~, u, l) (u-l)*cos(x)^2+l;
                obj.fvalinv = @(p, u, l) acos(sqrt((p-l)/(u-l)));
                obj.fgrad = @(x, u, l) -sin(2*x)*(u - l);
            else
                obj.fval = @(x, ~, ~, ~) x;
                obj.fvalinv = @(p, ~, ~) p;
                obj.fgrad = @(~, ~, ~) 1;
            end
            
            obj.ubound = @(~, u, ~) u;
            obj.uboundgrad = @(~, ~, N) zeros(1,N);
            obj.lbound = @(~, ~, l) l;
            obj.lboundgrad = @(~, ~, N) zeros(1,N);
        end % of constructor
        
        removeConstraint(obj, name);
        set(obj, varargin);
        n = getVarNum(obj, filter);
        p = map(obj, x);
        x = invmap(obj, p);
        names = getName(obj);
        jac = jacobian(obj, x);
        constraintList = getConstraints(obj);
        obj = vertcat(varargin);
        setName(obj, str);
        horzcat(varargin);
        %\\
    end % of methods (public)
end