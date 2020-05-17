classdef linker < handle
    % LINKER 
    %
    %   A LINKER object allows mapping an unconstrained optimisation 
    %   variable x to constrained model parameters p to be used in 
    %   conjuction with COMPARTMENT classes. 
    %
    %
    %   properties:
    %       name              - method for linking function
    %       upperBound        - upper bound for parameter p
    %       lowerBound        - lower bound for parameter p
    %       x                 - memory to store current x
    %       p                 - memory to store current p
    %
    %   methods (public):
    %       link              - forward transform p = F(x)  
    %       invLink           - backward transform x = F^(-1)(p)
    %       grad              - gradient of forward transform wrt to x
    %
    % Mohsen Farzi
    % Email: m.farzi@leeds.ac.uk
    
    properties (SetAccess = 'private')
        type = 'linear';      % Default is bounded cosine transform
        name = [];            % parameter name when printing messages
%         constrained = 'none'; % five types are constrained are possible in
%                               % LINKER class:
%                               % equality: Pi = Pj
%                               % inequality: Pi >= Pj
%                               % cumsum2: Pi + Pj = 1
%                               % cumsum3: Pi + Pj + Pk = 1
%                               % fixed: Pi = constant
        upperBound =  inf;
        lowerBound = -inf; 
        initRange = [0, 1];   % range for initialising the parameters
        compOrder = 1;        % computation order of linking parameters;
                              % integer number from 1 to N
    end
    
    properties (Access = 'private')                              
        constant = false;     % true if the parameter is constant at a 
                              % fixed value 
                           
        dummy = false;        % true if it is a dummy parameter that can be 
                              % experessed in terms of other model 
                              % parameters
                              %
                              % NOTE 1: A DUMMY LINK CANNOT SERVE AS A REF
                              % LINK FOR ANOTHER EQUALITY CONSTRAINT.
                              %
                              % NOTE 2: A DUMMY LINK CANNOT BE USED IN AN
                              % INEQUALITY CONSTRAINT AT ALL. USE THE
                              % CORRESPONDING REF LINK INSTEAD.
        
        linked = 0;           % The parameter is linked to another one
                              % 0 -> free parameter; no linking
                              % 1 -> being smaller than another parameter
                              % 2 -> cumsum 3b
        
        crossLinkNo = [];     % In case a parameer is dependent on others, 
                              % corresponding ref links are stored here.
                              % constraint.
                              
        storedParamVal = 0;   % parametert p = link(x). This parameter 
                              % serves as memory to track constant 
                              % parameters in model.    
                              
        fval = [];            % function handle to the mapping function for
                              % equality constraint
                        
        fvalinv = [];         % function handle to the inverse mapping 
                              % function                      
        
        gval = [];            % function handle to the gradient function of
                              % the mapping for the equality constraint
                              
        symbolicFun = [];       % symbolic operators for use in printing
                              % constraints                      
    end
    
    methods 
        function obj = linker(varargin)
            %LINKER Construct Function.
            %
            %   linker(..., <Parameter>, <Value>, ...) allow
            %   passing model parameters as parameter-value pairs.
            %           
            
            obj.set(varargin{:});
        end % of constructor
        
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
                    case 'squared'
                        obj.fval = @(x) x^2;
                        obj.fvalinv = @(p) sqrt(p);
                        obj.lowerBound = 0;
                        obj.upperBound = inf;
                        
                    case 'linear'
                        obj.fval = @(x) x;
                        obj.fvalinv = @(p) p;
                        obj.lowerBound = -inf;
                        obj.upperBound = inf;
                        
                    case 'cos'
                        obj.fval = @(x,ub,lb) cos(x)^2*(ub-lb)+ lb;
                        obj.fvalinv = @(p,ub,lb) acos(sqrt((p - lb)/(ub - lb)));
                end
            end % of if nargin==0
        end % of set
        
        function p = link(obj, x)            
            % check object format
            objIsNotColumnVector = not(iscolumn(obj));
            if objIsNotColumnVector
                error('Linker object must be a column vector or a scalar.');
            end
                                  
            % if required, extend input optimisation paramters with dummy
            % zero values
            objLen = size(obj, 1);
            if length(x) < objLen
                xPrime = zeros(objLen, 1);    
                isNotConstantOrDummy = not(obj.getConstant | ...
                                                  obj.getDummy);
                xPrime(isNotConstantOrDummy) = x;
            else
                xPrime = x;
            end
            
            [~,thisCompOrder] = sort(obj.getCompOrder());
            
            % initialise output parameter vector
            isNotComputed = true(objLen, 1);
            p = zeros(objLen, 1);
            
            for n = thisCompOrder'
                p(n) = obj(n).scalarLink(xPrime(n), p);
                isNotComputed(n) = false;
                
%                 if strcmp(obj(n).constrained, 'cumsum3a')
%                     linksNo = [n, obj(n).crossLinkNo];
%                     R = obj(n).c - sum(obj(linksNo).getLowerBound);
%                     p(n) = R * cos(xPrime(n))^2 + obj(n).lowerBound;
%                     isNotComputed(n) = false;
%                 end
%                 
%                 if strcmp(obj(n).constrained, 'cumsum3b')
%                     linksNo = [obj(n).crossLinkNo(1), n, obj(n).crossLinkNo(2)];
%                     R = obj(n).c - sum(obj(linksNo).getLowerBound);
%                     p(n) = R * sin(xPrime(linksNo(1)))^2 * cos(xPrime(n))^2 + obj(n).lowerBound; 
%                     isNotComputed(n) = false;
%                 end 
            end
            
            if any(isNotComputed)
                idx = find(isNotComputed);
                errMsg = [];
                for i = idx'
                    thisMsg = sprintf(['Error in link method for ', ...
                                       'links %d.\n'], i);
                     errMsg = strcat(errMsg, thisMsg);              
                end
                error(errMsg);
            end
        end
        
        function x = invLink(obj, p)
            if any(isnan(p))
                idx = find(isnan(p));
                wrnMsg = [];
                for i = idx
                    thisMsg = sprintf(['Input parameter has NaN ',...
                                       'values in links %d'], i);
                    wrnMsg = strcat(wrnMsg, thisMsg);
                end
                warning(wrnMsg);
                
                % replace nan values with random values
                p = obj.init();
            end
            dim = size(obj);
            if size(obj, 2)>1 || length(dim)>2
                error('Object should be a column vector or a scalar.');
            end

            % assert that obj dimension is the same as input p
            matchedDim = all(size(obj) == size(p));
            if ~matchedDim
                error(['The array object size should match the ', ...
                       'input array P.']);
            end

            % assert input p is in the range [lowerBound, upperBound]
            if any(p > obj.getUpperBound)
                idx = find(p > obj.getUpperBound);
                errMsg = [];
                for i = idx
                    if isempty(obj(i).name)
                        errMsg = strcat(errMsg, ...
                            sprintf(['Input parameter (%1.2e) is higher than ', ...
                            'the upper limit (%1.2e) for links %d.\n'], p(i), obj(i).upperBound, i));
                    else
                        errMsg = strcat(errMsg, ...
                            sprintf(['Input parameter (%1.2e) is higher than ', ...
                            'the upper limit (%1.2e) for "%s".\n'], p(i), obj(i).upperBound, obj(i).name));
                    end
                    % set p equal to the upper bound
                    p(i) = obj(i).upperBound;
                end
                warning(errMsg);
                %\\
            elseif any(p < obj.getLowerBound)
                idx = find(p < obj.getLowerBound);
                errMsg = [];
                for i = idx'
                    if isempty(obj(i).name)
                        errMsg = strcat(errMsg, ...
                            sprintf(['Input parameter (%1.2e) is lower than ', ...
                            'the lower limit (%1.2e) for links %d.\n'], p(i), obj(i).lowerBound, i));
                    else
                        errMsg = strcat(errMsg, ...
                            sprintf(['Input parameter (%1.2e) is lower than ', ...
                            'the lower limit (%1.2e) for "%s".\n'], p(i), obj(i).lowerBound, obj(i).name));
                    end
                    % set p equal to the lower bound
                    p(i) = obj(i).lowerBound;
                end
                warning(errMsg);
            end

            % compute the inverse link
            x = zeros(dim);
                                       
            % step3: check if alias links have the same values
%             isAlias = obj.getAlias();
%             crossLinks = arrayfun(@(thisObj) thisObj.crossLinkNo, obj(isAlias));
%             if any(p(isAlias) ~= p(crossLinks))
%                 idx = find(p(isAlias) ~= p(crossLinks));
%                 aliasLinkNo = find(isAlias);
%                 errMsg = [];
%                 for i = idx'
%                     ii1 = aliasLinkNo(i);
%                     ii2 = obj(ii1).crossLinkNo;
%                     thisMsg = sprintf(['Input parameters %d (%1.4f) does not ', ...
%                                        'equal links %d (%1.4f).\n'], ...
%                                        ii1, p(ii1), ii2, p(ii2));
%                     errMsg = strcat(errMsg, thisMsg);
%                 end
%                 error(errMsg);
%             end
            
            [~,thisCompOrder] = sort(obj.getCompOrder());
            % step4: compute values for inequality constraint
            for n = thisCompOrder'
                % values for constant links will be removed later
                x(n) = obj(n).scalarInvLink(p(n), p);
                
%                 % values for ineqaulity
%                 if strcmp(obj(n).constrained, 'inequality')
%                     idxGreater = obj(n).crossLinkNo;
%                     errMsg = sprintf(['Inequality is not valid ', ...
%                                       'for the constrained link %d'], n);
%                     assert(all(p(n) <= p(idxGreater)), errMsg);
%                     nom = p(n) - obj(n).lowerBound;
%                     denom = min(p(idxGreater(1)), obj(n).upperBound) - obj(n).lowerBound;
%                     x(n) = acos(sqrt(nom/denom));
%                 end
%                 
%                 if strcmp(obj(n).constrained, 'cumsum2a')
%                     linkNo = obj(n).crossLinkNo;
%                     denom = obj(n).c - obj(n).lowerBound - obj(linkNo).lowerBound;
%                     nom = p(n)-obj(n).lowerBound;
%                     x(n) = acos(sqrt(nom/denom));
%                 end
%                 
%                 if strcmp(obj(n).constrained, 'cumsum3a')
%                     linksNo = [n, obj(n).crossLinkNo];
%                     denom = obj(n).c - sum(obj(linksNo).getLowerBound);
%                     nom = p(n)-obj(n).lowerBound;
%                     x(n) = acos(sqrt(nom/denom));
%                 end
%                 
%                 if strcmp(obj(n).constrained, 'cumsum3b')
%                     linksNo = [obj(n).crossLinkNo(1), n, obj(n).crossLinkNo(2)];
%                     R = obj(n).c - sum(obj(linksNo).getLowerBound);
%                     denom = R - p(linksNo(1)) + obj(linksNo(1)).lowerBound;
%                     nom = p(n)-obj(n).lowerBound;
%                     x(n) = acos(sqrt(nom/denom));
%                 end
            end
            
            % remove dummy zeros values
            isConstantOrDummy = obj.getConstant | obj.getDummy;
            x(isConstantOrDummy) = [];
        end
        
        function jac = jacobian(obj, x)
            % grad(obj, x) return the gradient of link function wrt
            % optimisation parameter "x".
            dim = size(obj);
            if dim(2) > 1 || length(dim)>2
                error('Method "jacobian" only works for column arrays.');
            end
            
            thisP = obj.link(x);
            N = dim(1);
            jac = eye(N);
            
            % check if dummy or constant linkes exist
            if length(x) < N
                xPrime = zeros(N, 1);         % extended x vector with
                                              % dummy zero values
                isNotAliasOrConstantOrCumsum = not(obj.getConstant | obj.getDummy);
                xPrime(isNotAliasOrConstantOrCumsum) = x;
                x = xPrime;
            end
            
            [~,thisCompOrder] = sort(obj.getCompOrder());
            
            for n = thisCompOrder'
                jac(n,:) = obj(n).grad(x(n), n, thisP, jac);
            end
                
            % remove dummy zeros values
            isAliasOrConstantOrCumsum = obj.getConstant | obj.getDummy;
            jac(:,isAliasOrConstantOrCumsum) = [];
        end % of jacobian
        
            
        function thisCompOrder = getCompOrder(obj)
            thisCompOrder = arrayfun(@(thisObj) thisObj.compOrder, obj);
        end
        
        function setCompOrder(obj, thisCompOrder)
            for i=1:size(obj, 1)
                obj(i).compOrder = thisCompOrder(i);
            end
        end
    
        function v = getUpperBound(obj)
            v = arrayfun(@(thisObj) thisObj.upperBound, obj);
        end
        
        function v = getLowerBound(obj)
            v = arrayfun(@(thisObj) thisObj.lowerBound, obj);
        end
        
        function v = getConstant(obj)
            v = arrayfun(@(thisObj) thisObj.constant, obj);
        end
        
        function v = getDummy(obj)
            v = arrayfun(@(thisObj) thisObj.dummy, obj);
        end
        
        function constraintList = getConstraints(obj)
            %
             linkTypes = arrayfun(@(c) c.type, obj, 'UniformOutput', false);
             constraintList{1,1} = ['link type: ', strjoin(linkTypes, ',')];
            for idx = 1:length(obj)
                if ~isinf(obj(idx).upperBound) && ~obj(idx).constant && ~obj(idx).dummy
                    constraintList = [constraintList; sprintf('%s<=%1.2e', ...
                           obj(idx).name, obj(idx).upperBound)];
                end

                if ~isinf(obj(idx).lowerBound) && ~obj(idx).constant && ~obj(idx).dummy
                    constraintList = [constraintList; sprintf('%s>=%1.2e', ...
                        obj(idx).name, obj(idx).lowerBound)];
                end

                if obj(idx).constant
                    constraintList = [constraintList; sprintf('%s=%1.2e', ...
                                obj(idx).name, obj(idx).storedParamVal)];
                end

                if obj(idx).dummy
                    refLinks = obj(idx).crossLinkNo;
                    names = {obj(idx).name};
                    for i = 1:length(refLinks)
                        names = [names, {obj(refLinks(i)).name}];
                    end
                    constraintList = [constraintList; obj(idx).symbolicFun(names)];
                end

                if obj(idx).linked == 1
                    idx1 = idx;
                    idx2 = obj(idx).crossLinkNo;
                    constraintList = [constraintList; sprintf('%s>=%s', obj(idx2).name, obj(idx1).name)];
                end
            end
        end
        
        function setName(obj, str)
            for idx = 1:length(obj)
                obj(idx).name = str{idx};
            end
        end
        
        function p = init(obj)
            dim = size(obj);
            if dim(2)>1 || length(dim) > 2
                error(['Linker object should be either ', ...
                       'scalar or a column vector']);
            end
            
            p = zeros(dim);
            x = zeros(dim);
            for n = 1:dim(1)
                U = min(obj(n).upperBound, obj(n).initRange(2));
                L = max(obj(n).lowerBound, obj(n).initRange(1));
                p(n) = rand(1)*(U - L) + L;
                x(n) = obj(n).scalarInvLink(p(n));
            end
            
            % remove dummy or constant values
            isConstantOrDummy = obj.getConstant | obj.getDummy;
            x(isConstantOrDummy) = [];
            
            p = obj.link(x); 
        end % of init
        
        function obj = vertcat(varargin)
            argsNum = length(varargin);
            objectsNum = zeros(argsNum, 1); 
            
            % check if all inputs are of type linker and column-wise
            for nArgs = 1:argsNum
                thisObj = varargin{nArgs};
                if ~isa(thisObj, 'linker')
                    error(['Conversion from %s to linker', ...
                           ' is not possible.'], class(thisObj));
                end
                
                if size(thisObj,2) ~= 1
                    error('Linker class must be a column-wise object.');
                end
                
                objectsNum(nArgs) = size(thisObj, 1);
            end
            
            % 
            obj = linker.empty(sum(objectsNum), 0);
            i1 = 1;
            for nArgs = 1:argsNum
                i2 = i1 + objectsNum(nArgs) - 1;
                obj(i1:i2,1) = varargin{nArgs}.refreshLinks(i1-1);
                i1 = i2+1;
            end
            obj.sortCompOrder();
        end
        
        function horzcat(varargin)
            error(['Horizental concatation is not ', ...
                   'allowed for class linker']);
        end
        
%         function varargout = subsref(obj, s)
%             switch s(1).type
%               case '.'
%                  if length(s) == 1 && ~ismember(s.subs, methods(linker))
%                     % Implement obj.PropertyName
%                     ...
%                  elseif length(s) == 2 && strcmp(s(2).type,'()')
%                     % Implement obj.PropertyName(indices)
%                     ...
%                  else
%                     [varargout{1:nargout}] = builtin('subsref',obj,s);
%                  end
%               case '()'
%                  if length(s) == 1
%                     % Implement obj(indices)
%                     ...
%                  elseif length(s) == 2 && strcmp(s(2).type,'.')
%                     % Implement obj(ind).PropertyName
%                     ...
%                  elseif length(s) == 3 && strcmp(s(2).type,'.') && strcmp(s(3).type,'()')
%                     % Implement obj(indices).PropertyName(indices)
%                     ...
%                  else
%                     % Use built-in for any other expression
%                     [varargout{1:nargout}] = builtin('subsref',obj,s);
%                  end
%               case '{}'
%                  if length(s) == 1
%                     % Implement obj{indices}
%                     ...
%                  elseif length(s) == 2 && strcmp(s(2).type,'.')
%                     % Implement obj{indices}.PropertyName
%                     ...
%                  else
%                     % Use built-in for any other expression
%                     [varargout{1:nargout}] = builtin('subsref',obj,s);
%                  end
%               otherwise
%                  error('Not a valid indexing expression')
%            end
%         end % of subsref
        %\\
    end % of methods (public)
    
    methods (Access = 'private')
%         function setScalarX(obj, x)
%             obj.x = x;
%         end
        
%         function thisX = getScalarX(obj)
%             if strcmp(obj.type, 'dummy')
%                 thisX = obj.aliasLink.getX;
%             else
%                 thisX = obj.x;
%             end
%         end

        function obj = refreshLinks(obj, offset)
            for n = 1:size(obj,1)
                obj(n).compOrder = obj(n).compOrder + offset;
                if ~isempty(obj(n).crossLinkNo)
                    obj(n).crossLinkNo = obj(n).crossLinkNo + offset;
                end
            end
        end 
        
        function sortCompOrder(obj)
            % sort the computation order following concatanation
            thisCompOrder = obj.getCompOrder();
            N = size(obj, 1);
            
            % sort constant parameters at the top of the list, i.e. 
            % constant parameters should be addressed first
            isConstant = obj.getConstant();
            n = find(isConstant);
            for i = 1:length(n)
                idx = thisCompOrder < thisCompOrder(n(i));
                thisCompOrder(idx) = thisCompOrder(idx) + 1;
                thisCompOrder(n(i)) = 1;
            end
            
            % sort dummy parameters at the end of the list; i.e. dummy
            % values should be computed once all other params are computed.
            isDummy = obj.getDummy();
            n = find(isDummy);
            for i = 1:length(n)
                idx = thisCompOrder > thisCompOrder(n(i));
                thisCompOrder(idx) = thisCompOrder(idx) - 1;
                thisCompOrder(n(i)) = N;
            end
            
            % update compuation order
            obj.setCompOrder(thisCompOrder);
        end
        
        function thisP = scalarLink(obj, thisX, p)
            % fun(obj, x) define model parameter "p" in terms of the 
            % optimisation parameter "x". 
            %
            % see also: invFun, grad
            %
            if nargin == 1
                thisX = 0;
                p = 0;
            elseif nargin == 2
                p = 0;
            end
            
            switch obj.type
                case 'cos'
                    lBound = obj.lowerBound;
                    if obj.linked == 0
                        uBound = obj.upperBound;
                    elseif obj.linked == 1
                        uBound = p(obj.crossLinkNo);
                    elseif obj.linked == 2
                        uBound = obj.upperBound - p(obj.crossLinkNo);
                    else
                        error('Unrecognised value %d for linked property', obj.linked);
                    end
                    thisP = cos(thisX)^2*(uBound-lBound)+lBound;
 
                case 'squared'
                    thisP = thisX^2;
                        
                case 'linear'
                    thisP = thisX;
                
                case 'dummy'
                    % note for dummy varialbes, thisX is an alias name for
                    % thisP
                    thisP = obj.fval(p, obj.crossLinkNo);
                    
                case 'constant'
                    thisP = obj.storedParamVal;
                    
                 otherwise
                    error('The linking type %s is not recognised.', obj.type);
            end % switch obj.type
        end % of scalarLink method
        
        function thisX = scalarInvLink(obj, thisP, p)
            % invFun(obj, p) define the optimisation parameter "x" in 
            % terms of the model parameter "p". 
            %
            % see also: invFun, grad
            %
            switch obj.type
                case 'cos'
                    lBound = obj.lowerBound;
                    if obj.linked == 0 || nargin == 2
                        uBound = obj.upperBound;
                    elseif obj.linked == 1
                        uBound = p(obj.crossLinkNo);
                    elseif obj.linked == 2
                        uBound = obj.upperBound - p(obj.crossLinkNo);
                    else
                        error('Unrecognised value %d for linked property', obj.linked);
                    end
                    
                    thisX = acos(sqrt((thisP - lBound)/(uBound - lBound)));
                
                case 'squared'
                    thisX = sqrt(thisP);
                    
                case 'linear'
                    thisX = thisP;
                
                case 'dummy'
                    thisX = 0;
                    
                case 'constant'   
                    thisX = 0;
                    
                 otherwise
                    error('The linking type %s is not recognised.', obj.type);
            end % switch obj.type
        end
        
        function g = grad(obj, x, n, p, thisJac)
            nObj = length(p);
            
            switch obj.type
                case 'cos'
                    lBound = obj.lowerBound;
                    if obj.linked == 0
                        uBound = obj.upperBound;
                        g = zeros(1, nObj);
                        g(n) = -sin(2*x)*(uBound - lBound);
                        
                    elseif obj.linked == 1
                        g = cos(x)^2*thisJac(obj.crossLinkNo, :);
                        uBound = p(obj.crossLinkNo);
                        g(n) = -sin(2*x)*(uBound - lBound);
                        
                    elseif obj.linked == 2
                        g = -cos(x)^2*thisJac(obj.crossLinkNo, :);
                        uBound = obj.upperBound - p(obj.crossLinkNo);
                        g(n) = -sin(2*x)*(uBound - lBound);
                    else
                        error('Unrecognised value %d for linked property', obj.linked);
                    end
                    
                    
                case 'squared'
                    g = zeros(1, nObj);
                    g(n) = 2*x;
                    
                case 'linear'
                    g = zeros(1, nObj);
                    g(n) = 1;
                    
                case 'dummy'
                    g = obj.gval(thisJac, p, obj.crossLinkNo);
                    
                case 'constant'    
                    g = zeros(1, nObj);
                otherwise
                    error('The linking type %s is not recognised.', obj.type);
           end
        end % of grad
        
        function isValidInputForCumsum2(obj, link1No, link2No, c)
            % check the validity of input paramters
            N = length(obj);
            
            % check if link1 is a valid link number
                assert(link1No>0 && link1No<N+1, ...
                       sprintf(['Input link should be an interger ', ...
                                'between 1 and %d.'], length(obj)));
                            
                % check if link2 is a valid link number
                assert(link2No>0 && link2No<N+1, ...
                       sprintf(['Input link should be an interger ', ...
                                'between 1 and %d.'], length(obj)));
                
                % check if c is between 0 and 1
                assert(c>=0 && c<=1, ...
                       'Sum of squared links should be between 0 and 1');
                   
                link1TypeIsNotCos = ~strcmp(obj(link1No).type, 'cos');
                link2TypeIsNotCos = ~strcmp(obj(link2No).type, 'cos');
                if link1TypeIsNotCos || link2TypeIsNotCos
                    error(['For cumsum2 constraints, both ', ...
                           'arguments should be of tyep "cos".']);
                end
                
                if obj(link1No).dummy || obj(link2No).dummy
                    error('Input links cannot be dummy to other links.');
                end
                
                if obj(link1No).linked>0|| obj(link2No).linked>0
                    error(['Inout links have been used ', ...
                           'in an inequality constraint.']);
                end
        end
        
        function isValidInputForCumsum3(obj, link1No, link2No, link3No, c)
            % check the validity of input paramters
            N = length(obj);
            
            % check if link1 is a valid link number
                assert(link1No>0 && link1No<N+1, ...
                       sprintf(['Input link should be an interger ', ...
                                'between 1 and %d.'], length(obj)));
                            
                % check if link2 is a valid link number
                assert(link2No>0 && link2No<N+1, ...
                       sprintf(['Input link should be an interger ', ...
                                'between 1 and %d.'], length(obj)));
                
                % check if link3 is a valid link number
                assert(link3No>0 && link3No<N+1, ...
                       sprintf(['Input link should be an interger ', ...
                                'between 1 and %d.'], length(obj)));
                
                % check if link numbers are different
                assert(link1No ~= link2No, ['The first and the second ',...
                      'links are the same in the cumsum3 constraint.']);
                  
                assert(link1No ~= link3No, ['The first and the third ',...
                      'links are the same in the cumsum3 constraint.']);  
                  
                assert(link2No ~= link3No, ['The second and the third ',...
                      'links are the same in the cumsum3 constraint.']);  
                  
                % check if c is between 0 and 1
                assert(c>=0 && c<=1, ['Constant value should be ',...
                      'between 0 and 1 in cumsum3 constraint.']);
                   
                link1TypeIsNotCos = ~strcmp(obj(link1No).type, 'cos');
                link2TypeIsNotCos = ~strcmp(obj(link2No).type, 'cos');
                if link1TypeIsNotCos || link2TypeIsNotCos
                    error(['In cumsum3 constraints, all ', ...
                           'links should be of tyep "cos".']);
                end
                
                if obj(link1No).dummy || obj(link2No).dummy || obj(link3No).dummy
                    error('Input links cannot be dummy to other links.');
                end
                
                if obj(link1No).linked==1 || ...
                   obj(link2No).linked==1 || ...
                   obj(link3No).linked==1
                    error(['Inout links have been used ', ...
                           'in an inequality constraint.']);
                end
        end % of isValidInputForCumsum3
        
    end % of methods (private)
    
    methods (Static)
        function state = isValidLinker(v)
            state = isa(v, 'boundedCosLinker');
        end
        
    end
end

function state = isValidLinkingType(v)
            state = ismember(v, {'cos', 'squared', 'linear'});
end