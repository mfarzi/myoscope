function img = ifft3c(K)
% inverse fft transform for each diffusion weighted MR signal
img = zeros(size(K));

for i = 1:size(K,4) 
    thisKspace = K(:,:,:,i);
    img(:,:,:,i) = sqrt(length(thisKspace(:)))*ifftn(ifftshift(thisKspace));
end