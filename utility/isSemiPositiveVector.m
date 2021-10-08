function bool = isSemiPositiveVector(v)
    bool = isempty(v) || (isnumeric(v) && all(v>=0) && isvector(v));
end