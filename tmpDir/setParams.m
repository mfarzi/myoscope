function setParams(obj, varargin)
% setParams is a method for class COMPARTMENT
%
%   setParams(obj, p) update the class property "params" with the input
%   vector p.
%
%   setParams(obj, <paramName>, <paramValue>) update specific parameters
%   using params-value pairs. See getParamsName for a full list of model 
%   parameters.
%

% check the number of input arguments
assert(nargin>1, 'MATLAB:compartment:setParams', ...
    'Too few input parameters.');

if nargin==2
    % setParams(obj, p)
    p = varargin{1};
else
    % setParams(obj, <paramName>, <paramValue>)
    p = readPairs(obj, varargin{:});
end

% obj.validateParams(p);
obj.params = p;
end

function p = readPairs(obj, varargin)
% helper function to parse input params-value pairs

% default parameters
p = obj.params;

% get model parameter names
paramsList = obj.getParamsName;

% update specific parameters 
for i=1:2:nargin-1
    thisParam = varargin{i};
    assert(isa(thisParam, 'char'), 'MATLAB:compartment:setParams',...
        'Input parameter name must be of type char');
    idx = strcmp(paramsList, thisParam);
    if any(idx)
        v = varargin{i+1};
        validateattributes(v, {'numeric'}, {'nonnegative','scalar'},...
            'compartment.setParams', paramsList{i});
        p(idx) = v;
    else
        error('MATLAB:compartment:setParams',...
            '"%s" is not a valid parameter name.', thisParam);
    end
end
end