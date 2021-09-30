function bool = isDirectionMatrix(x)
    bool = ~isempty(x) && isnumeric(x) && ismatrix(x) && size(x, 2)==3;
end