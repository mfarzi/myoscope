function addDiffusionTimes(obj, dt)
    % addDiffusionTimes is a public method for class SCHEME
    %
    % addDiffusionTimes set discrete diffusion times.
    assert(isPositiveVector(dt), 'MATLAB:scheme:addDiffusionTimes',...
        'Input argument must be a positive vector.');
    
    if isrow(dt)
        dt = dt';
    end
    
    obj.diffusionTimes = union(obj.diffusionTimes, dt, 'stable');
end
    
    