function b = computeBvalue(gmag, dt, delta)
    % estiamte b-value
    b = (dt-delta/3).*((delta .*gmag)*math.GAMMA).^2;
end