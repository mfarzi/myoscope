classdef optimizer < handle
    % OPTIMISER 
    %
    %   An OPTIMISER object encapsulates various numerical techniques for 
    %   finding a local minimum by minimizing a scalar cost function f(x) 
    %   with respect to the vector of parameters x. 
    %   
    %   NOTE: the objective function and its derivatives must be 
    %   implemented in the subclasses.
    %
    %   For mathematical background see
    %       Jorge Nocedal and Stephen J. Wright, "Numerical Optimization",
    %       Springer Series in Operation Research and Financial Engineering
    %       , 2006.
    %
    %   optimizer methods:
    %       run           - start the optimisation procedure from x0
    %       setInitPoint  - set the initial point x0 to start the
    %                       optimisation
    %       setOptParams  - set the optimisation parameters
    %
    %   See also: compartment
    
    % Mohsen Farzi
    % Email: m.farzi@leeds.ac.uk
    
    properties (SetAccess = 'private')
        % optimisation algorithms: 
        %       'levenberg-marquardt' (default)
        %       'conjugate-gradient'
        algorithm = 'levenberg-marquardt';

        % Level of display
        % false: displays no output (default).
        % true : displays output at each iteration
        prettyPrint = false;
        
        % stop tollerance for variation in cost function
        fevalTol = 1e-6;
        
        % stop tollerance for variation in function gradient
        fgradTol = 1e-6;
        
        % maximum number of iterations for the main algorithm
        maxIter = 500;                  
    end
    
    properties (Hidden)
        % wolf line search: parameter used in bracketing phase to estimate
        % a_j and b_j at [0,alpha_max]
        wolfLineSearch_alphaMax = 100;                   

        % Sufficient Decrease Condition (WC1) 
        % typical values: 0.01 or 1e-4
        % f(x_k + alpha*P_k) <= f(x_k) + c1*alpha*P_k'*gf(x_k) (WC1)
        wolfLineSearch_c1 = 1e-4;       

        % parameter for WC2 with typical value
        % of 0.8 or 0.9
        % |P_k'.gf(x_k + alpha*P_k)| <= c2|P_k'*gf(x_k)|       (WC2) 
        wolfLineSearch_c2 = 0.8;       
    end
    
    properties(Access = 'private')
        feval = [];                 % handle to cost function
        fjac = [];                  % handle to gradient
    end
    
    methods 
        function obj = optimizer(varargin)
            %OPTIMISER Construct Function.
            %   optimizer() constructs an empty object with default set of
            %   parameters. 
            %   NOTE: use setfeval(@fhandle) or setfjac(@fhandle) to sef 
            %   function handles later.
            %
            %   optimizer(f, gf) constructs an object with function handles
            %   f and gf to with default sets of parameters.
            %
            %   optimizer(..., <Parameter>, <Value>, ...) allow
            %   passing more inforamtion as parameter-value pairs. 
            % 
            if nargin == 0
                % >> obj = optimizer();
                % Do Nothing!
            else
                % update the rest of parameters using the method "set"
                obj.set(varargin{:});
            end
        end % of construct function
        
        function set(obj, varargin)
            if nargin == 1
                % >> obj = obj.set()
                % Don Nothing!
            else
                p = inputParser;
                p.CaseSensitive = false;
                
                p.addParameter('algorithm' , obj.algorithm, @(v) ischar(v));
                p.addParameter('prettyPrint' , obj.prettyPrint, @(v) islogical(v) && isscalar(v));
                p.addParameter('maxIter'     , obj.maxIter , @(v) isnumeric(v) && isscalar(v));
                p.addParameter('fevalTol'         , obj.fevalTol , @(v) isnumeric(v) && isscalar(v));
                p.addParameter('fgradTol'         , obj.fgradTol , @(v) isnumeric(v) && isscalar(v));
                p.addParameter('c1'         , obj.wolfLineSearch_c1 , @(v) isnumeric(v) && isscalar(v) && v<1 && v>0);
                p.addParameter('c2'         , obj.wolfLineSearch_c2 , @(v) isnumeric(v) && isscalar(v) && v<1 && v>0);
                p.addParameter('feval'      , obj.feval , @(v) isa(v, 'function_hadle'));
                p.addParameter('fjac'      , obj.fjac , @(v) isa(v, 'function_hadle'));

                p.parse(varargin{:});
                
                obj.algorithm = p.Results.algorithm;
                obj.prettyPrint = p.Results.prettyPrint;
                obj.maxIter = p.Results.maxIter;
                obj.fevalTol = p.Results.fevalTol;
                obj.fgradTol = p.Results.fgradTol;
                obj.wolfLineSearch_c1 = p.Results.c1;
                obj.wolfLineSearch_c2 = p.Results.c2;
                obj.feval = p.Results.feval;
                obj.fjac = p.Results.fjac;
            end
        end
        
        function setfeval(obj, fhandle)
            if isa(fhandle, 'function_handle')
                obj.feval = fhandle;
            else
                error('Input should be a function handle.\n');
            end
        end
        
        function setfjac(obj, fhandle)
            if isa(fhandle, 'function_handle')
                obj.fjac = fhandle;
            else
                error('Input should be a function handle.\n');
            end
        end
        
        function [x, cost, exitFlag, output] = run(obj, x0, varargin)
            if ~isa(obj.feval, 'function_handle')
                error('feval is not set properly!.\n');
            end
            
            if ~isa(obj.fjac, 'function_handle')
                error('fjac is not set properly!\n');
            end
            
            if strcmp(obj.algorithm, 'conjugate-gradient')
                options = obj.getOptions();
                fval = @(x) sum(obj.feval(x, varargin{:}).^2);
                gval = @(x) 2*obj.fjac(x, varargin{:})'*obj.feval(x, varargin{:});
                [x, cost, exitFlag, output] = conjugateGradient(fval, gval, x0, options);
            elseif strcmp(obj.algorithm, 'levenberg-marquardt')
                options = obj.getOptions();
                objective = {@(x) obj.feval(x, varargin{:}), @(x) obj.fjac(x, varargin{:})};
                [x, cost, ~, exitFlag, output] = lsqnonlin(objective, x0, [], [], options);
            elseif strcmp(obj.algorithm, 'interior-point')
                options = obj.getOptions();
                fval = @(x) sum(obj.feval(x, varargin{:}).^2);
                gval = @(x) 2*obj.fjac(x, varargin{:})'*obj.feval(x, varargin{:});
                A = [[1,  -1,  0, 0, 0, 0];
                     [ 0, 1,  -1, 0, 0, 0]];
                b = [-1e-10; -1e-10];
                %Aeq = [[1, -1,  0, 0, 0, 0];
                %       [1,  0, -1, 0, 0, 0]];
                %beq = [0; 0];   
                Aeq = [];
                beq = [];
                lb = [];
                ub = [];
                [x,cost,exitFlag,output] = fmincon({fval, gval},x0,A,b,Aeq,beq,lb,ub, [], options);
            end
        end % of run
    end %of methods
end % of the class optimizer