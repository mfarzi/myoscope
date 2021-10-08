function bool = isDirectionMatrix(x)
    bool = isnumeric(x) && ismatrix(x) && size(x, 2)==3;
end