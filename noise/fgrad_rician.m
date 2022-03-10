function g = fgrad_rician(x, data)
    % 
    % input parameters:
    %              sig: Synthesised signal from the compartment model
    %             data: Measured DW-MR signal
    %            sigma: Rician noise parameter
    assert(iscolumn(data) && all(data>=0),...
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
    
    I1_I0_z = zeros(size(z));
    I1_I0_z(idx) = besseli(1,z(idx))./besseli(0,z(idx));
    I1_I0_z(~idx) = 1 - exp((z(~idx)-700)*(-0.00143010502604)-7.24386992534278);
    %L = sum(sig.^2)/(2*sigma2)-sum(log(besseli(0, z)));
    gL_gA = (N*A)/sigma2 - sum(data/sigma2.*(I1_I0_z));
    gL_gSigma =-(N*A^2)/(sigma*sigma2)+2/sigma*sum(z.*I1_I0_z)+...
               2*N/sigma-(data'*data)/(sigma*sigma2); 
    g = zeros(2,1);
    g(1) = gL_gA*2*x(1);
    g(2) = gL_gSigma*2*x(2);
end
