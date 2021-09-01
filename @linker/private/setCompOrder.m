function setCompOrder(obj, thisCompOrder)
    for i=1:size(obj, 1)
        obj(i).compOrder = thisCompOrder(i);
    end
end