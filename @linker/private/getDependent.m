function idx = getDependent(obj, index)
    % getDependent is a private method for class linker
    %
    % getDependent(obj, index) return all indices required to comput
    % function value at input index
    idx = index;
    for i = index
        if isempty(obj(i).dependence)
            continue;
        else
            idx = [idx, obj.getDependent(obj(i).dependence)];
        end
    end
end