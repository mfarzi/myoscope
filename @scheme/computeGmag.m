function g = computeGmag(bval, dt, delta)
    % estiamte gmag
    c = (dt-delta/3).*(delta.^2)*math.GAMMA^2;
    g = sqrt(bval./c);
end