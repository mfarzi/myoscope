function L = ricianLogLikelihood(data, sig, sigma)
    % 
    % input parameters:
    %              sig: Synthesised signal from the compartment model
    %             data: Measured DW-MR signal
    %            sigma: Rician noise parameter
    assert(iscolumn(data) && all(data>=0),...
        'MATLAB:noiseModel:invalidInputArgs',...
        'First input argument must be a non-negative column vector.');
    N = length(data);
    
    assert(iscolumn(sig) && length(sig)==N && all(sig>=0),...
        'MATLAB:noiseModel:invalidInputArgs',...
        ['Second input argument must be a non-negative column vector',...
         ' of size %d.'], N);
    
    assert(isscalar(sigma)&&sigma>0,...
        'MATLAB:noiseModel:invalidInputArgs',...
        'Third input argument must be a positive scalar.');
    
    % rician noise pdf
    %P = @(m, A, sigma2) (m/sigma2).*exp(-(A.^2+m.^2)/(2*sigma2)).*besseli(0,(A.*m)/sigma2);
    sigma2 = sigma^2;
    z = (sig.*data)/sigma2;
    idx = z<700;
    logI0_z = zeros(size(z));
    logI0_z(idx) = log(besseli(0, z(idx)));
    logI0_z(~idx) = z(~idx)-0.5*log(2*pi*z(~idx));
    L = (sig'*sig)/(2*sigma2)-sum(logI0_z)+(data'*data)/(2*sigma2)-sum(log(data/sigma2));
end
