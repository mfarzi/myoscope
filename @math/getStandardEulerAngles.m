function getStandardEulerAngles(obj)
    % rotateAxis remove the ambiguity in estimation of parameters
    % theta and phi.
    %
    % the direction of v1 or -v1 is selected such that the angle
    % theta is always in range [0, pi/2].
    %
    % see also: getUnitFrame, getEulerAngles

    obj.theta = mod(obj.theta, 2*pi);
    obj.phi = mod(obj.phi, 2*pi);

    if obj.theta > pi
        obj.theta = 2*pi - obj.theta;
        obj.phi = mod(pi+obj.phi, 2*pi);
    end

    % make sure theta < pi/2
    if obj.theta>pi/2 
        obj.theta = pi - obj.theta;
        obj.phi = mod(pi + obj.phi, 2*pi);
    end

    obj.modelParams(4:5) = [obj.theta; obj.phi];
end % of rotateAxis