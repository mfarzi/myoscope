function g = normaliseDirections(ghat)
    % a static method for class scheme
    assert(isDirectionMatrix(ghat), 'MATLAB:scheme:addNominalGhat',...
        'Input argument must be a matrix with 3 columns.');
    
    ZERO_NORM_TOL = 0.001; 
    
    gNorm = sqrt(sum(ghat.^2,2));
    isZero = gNorm<ZERO_NORM_TOL;
    g = ghat./gNorm;
    g(isZero,1) = 0;
    g(isZero,2) = 0;
    g(isZero,3) = 0;
end