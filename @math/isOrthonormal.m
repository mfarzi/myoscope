function state = isOrthonormal(U)
    % ceck if 
    EPS = 1e-10;
    validateattributes(U, {'numeric'}, {'square'},...
        'isOrthonormal', 'U');
    n = size(U,1);
    I = eye(n);
    Ip = U'*U;
    
    state = all(abs(Ip(:)-I(:))<EPS);
end