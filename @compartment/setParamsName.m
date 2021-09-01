function setParamsName(obj,varargin)
    % setParamsName is a method for class MULTICOMPARTMENT
    %
    %   setParamsName(obj, [indices], names) set smantic names for the given
    %   parameter indices.
    %           Input Args:
    %                  obj: A MULTICOMPARTMENT object
    %              indices: An optional positional argument. A vector of 
    %                       integers from 1 to N, where N is the total
    %                       number of model parameters.
    %                       (Default = [1:7])
    %                names: A list of strings with the same length as
    %                       indices.
    
    [indices, names] = parseInputArguments(obj, varargin{:});
    obj.links(indices).setName(names);
end

function [indices, names] = parseInputArguments(obj, varargin)
    if nargin<2
        error('MATLAB:MULTICOMPARTMENT:setParamsName',...
              'Too few input arguments.');
    elseif nargin==2
        indices = 1:obj.nParams;
        names = varargin{1};
    elseif nargin==3
        indices = varargin{1};
        names = varargin{2};
    else
        error('MATLAB:MULTICOMPARTMENT:setParamsName',...
              'Too many input arguments.');
    end
    
    % validate input attributes
    isPositiveInteger = isnumeric(indices) && all(mod(indices,1)==0) &&...
                        all(indices>0);       
    assert(isPositiveInteger&&isvector(indices)&&all(indices<=obj.nParams),...
        'MATLAB:MULTICOMPARTMENT:setParamsName',...
        'Indices must be a vector of integers with values between 1 and %d.',...
        obj.nParams);
    if iscolumn(indices)
        indices = indices';
    end
    
    isaCellofStrings = iscellstr(names)&&length(names)==length(indices);
    isaScalarString = ischar(names)&&length(indices)==1;
    assert(isaCellofStrings || isaScalarString,...
            'MATLAB:MULTICOMPARTMENT:setParamsName',...
            'Names must be an array of strings with %d elements.',...
            length(indices));
    if isaScalarString
        names = {names};
    end
end