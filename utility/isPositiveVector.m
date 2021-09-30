function bool = isPositiveVector(v)
    bool = isnumeric(v) && all(v>0) && isvector(v);
end