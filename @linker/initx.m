function x = initx(obj)
    % initx is a public method to initialise x
    N = obj.getVarNum();
    x = zeros(N, 1);
    for n=1:N
        if obj(n).bounded
            x(n) = acos(sqrt(rand(1)));
        else
            x(n) = rand(1)*2*pi;
        end
    end
    
    % remove dummy or constant variables
    isDummy = strcmp({obj.type}', 'dummy');
    isConstant = strcmp({obj.type}', 'constant');
    x(isConstant|isDummy) = [];
    %x = rand(N,1)*pi/2; 
end