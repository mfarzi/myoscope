function options = getOptions(obj)
    %
    if strcmp(obj.algorithm, 'conjugate-gradient')
        options.alphaMax = obj.wolfLineSearch_alphaMax;
        options.alpha1 = 0;
        options.c1 = obj.wolfLineSearch_c1;
        options.c2 = obj.wolfLineSearch_c2;
        options.stopTol = obj.fevalTol;
        options.maxIter = obj.maxIter;
        options.prettyPrint = obj.prettyPrint;
    end
    
    if strcmp(obj.algorithm, 'levenberg-marquardt')
        if obj.prettyPrint
            display = 'iter';
        else
            display = 'off';
        end
        options = optimoptions('lsqnonlin', ...
                               'Algorithm','levenberg-marquardt', ...
                               'SpecifyObjectiveGradient',true  , ...
                               'MaxIterations', obj.maxIter, ...
                               'Display', display                  , ...
                               'FunctionTolerance', obj.fevalTol);
    end     
    
    if strcmp(obj.algorithm, 'interior-point')
        if obj.prettyPrint
            display = 'iter';
        else
            display = 'off';
        end
        options = optimoptions('fmincon', ...
                               'Algorithm','interior-point', ...
                               'SpecifyObjectiveGradient',true  , ...
                               'MaxIterations', obj.maxIter, ...
                               'Display', display                  , ...
                               'FunctionTolerance', obj.fevalTol, ...
                               'ConstraintTolerance', 1e-15);
    end   
end
    
    