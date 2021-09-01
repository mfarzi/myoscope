function obj = vertcat(varargin)
    % vertcat is a public method for class LINKER
    %
    % vertcat(obj1, obj2, ...) concatenate a new copy of input objects in a 
    % column format.
    
    newObj = cell(1, nargin);
    offset = 0;
    for n = 1:nargin
        % check input object
        validateattributes(varargin{n}, {'linker'}, {'column'});
        % copy it
        newObj{n} = copy(varargin{n});
        % shift link numbers
        newObj{n}.shiftCompOrder(offset);
        % update offset
        offset = offset + size(newObj{n},1);
    end
    
    % concatenate input objects
    obj = cat(1, newObj{:});  
    N = size(obj,1);
    
    % update index and varNum
    [obj.varNum] = deal(N);
    indices = num2cell(1:N);
    [obj.index] = indices{:};
    
    % update computation order
    obj.updateCompOrder();
end