function bool = isPositiveVector(v)
    bool = isempty(v) || (isnumeric(v) && all(v>0) && isvector(v));
end