 function set(obj, varargin)
    % obj.set() do nothing!
    %
    % obj.set(..., <Parameter>, <Value>, ...) 
    % allow passing class properties as parameter-value pairs. 
    %
    if nargin==1
        % >> obj = obj.set()
        % Do Nothing!
    else
        % Parse the parameter-value pairs
        parser = inputParser;
        parser.CaseSensitive = false;
        parser.addParameter('upperBound'  , obj.upperBound , ...
                            @(v) isnumeric(v) && isscalar(v));
        parser.addParameter('lowerBound'  , obj.lowerBound , ...
                            @(v) isnumeric(v) && isscalar(v)); 
        parser.addParameter('initRange'  , obj.initRange , ...
                            @(v) isnumeric(v) && all(size(v) == [1,2]));                 
        parser.addParameter('type'  , obj.type , ...
                            @(v) isValidLinkingType(v)); 
        parser.addParameter('name'  , obj.name , ...
                            @(v) isa(v, 'string') || isa(v, 'char'));                 

        parser.parse(varargin{:});              

        obj.upperBound   = parser.Results.upperBound; 
        obj.lowerBound   = parser.Results.lowerBound; 
        obj.initRange    = parser.Results.initRange;
        obj.type         = parser.Results.type; 
        obj.name         = parser.Results.name;


        switch obj.type
            case 'linear'
                obj.fval = @(x,~,~) x;
                obj.fvalinv = @(p,~,~) p;
                obj.lbound = [];
                obj.ubound = [];

            case 'cos'
                obj.fval = @(x,ub,lb) cos(x)^2*(ub-lb)+ lb;
                obj.fvalinv = @(p,ub,lb) acos(sqrt((p - lb)/(ub - lb)));
        end
    end % of if nargin==0
end % of set