function [gL_gSig, gL_gSigma] = ricianLogLikelihoodJacobian(data, sig, sigma)
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
    
    I1_I0_z = zeros(size(z));
    I1_I0_z(idx) = besseli(1,z(idx))./besseli(0,z(idx));
    I1_I0_z(~idx) = 1 - exp((z(~idx)-700)*(-0.00143010502604)-7.24386992534278);
    %L = sum(sig.^2)/(2*sigma2)-sum(log(besseli(0, z)));
    gL_gSig = sig/sigma2 - data/sigma2.*(I1_I0_z);
    gL_gSigma =-(sig'*sig)/(sigma*sigma2)+2/sigma*sum(z.*I1_I0_z)+...
               2*N/sigma-(data'*data)/(sigma*sigma2); 
end
