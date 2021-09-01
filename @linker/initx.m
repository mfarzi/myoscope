function x = initx(obj)
    % initx is a public method to initialise x
    N = obj.getVarNum('free');
    x = rand(N,1)*pi/2; 
end