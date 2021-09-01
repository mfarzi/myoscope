function state = isOrthonormal(U)
    % ceck if 
    validateattributes(U, {'numeric'}, {'square'},...
        'isOrthonormal', 'U');
    n = size(U,1);
    I = eye(n);
    Ip = U'*U;
    
    state = all(abs(Ip(:)-I(:))<eps);
end