function name = getName(obj)
name = arrayfun(@(c) c.name, obj, 'UniformOutput', false);
end