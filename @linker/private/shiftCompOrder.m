function obj = shiftCompOrder(obj, offset)
    % shiftCompOrder is a private method for class LINKER
    %
    % shiftCompOrder add offset value to link IDs
    objLen = size(obj, 1);
    for n = 1:objLen
        obj(n).compOrder = obj(n).compOrder + offset;
        if ~isempty(obj(n).dependence)
            obj(n).dependence = obj(n).dependence + offset;
        end
    end
end 