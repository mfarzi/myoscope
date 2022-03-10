function L = feval_rician(x, data)
    % 
    % input parameters:
    %              sig: Synthesised signal from the compartment model
    %             data: Measured DW-MR signal
    %            sigma: Rician noise parameter
    assert(iscolumn(data) && all(data>0),...
        'MATLAB:noiseModel:invalidInputArgs',...
        'First input argument must be a non-negative column vector.');
    N = length(data);
    
    A = x(1)^2+0.0001;
    sigma = x(2)^2+0.0001;
    
    % rician noise pdf
    %P = @(m, A, sigma2) (m/sigma2).*exp(-(A.^2+m.^2)/(2*sigma2)).*besseli(0,(A.*m)/sigma2);
    sigma2 = sigma^2;
    z = (A*data)/sigma2;
    idx = z<700;
    logI0_z = zeros(size(z));
    logI0_z(idx) = log(besseli(0, z(idx)));
    logI0_z(~idx) = z(~idx)-0.5*log(2*pi*z(~idx));
    L = (N*A^2)/(2*sigma2)-sum(logI0_z)+(data'*data)/(2*sigma2)-sum(log(data/sigma2));
end