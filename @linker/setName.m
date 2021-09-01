function setName(obj, str)
    for idx = 1:length(obj)
        obj(idx).name = str{idx};
    end
end