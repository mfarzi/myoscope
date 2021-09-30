function bool = isSemiPositiveVector(v)
    bool = isnumeric(v) && all(v>=0) && isvector(v);
end