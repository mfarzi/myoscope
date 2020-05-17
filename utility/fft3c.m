function K = fft3c(img)
% compute 3D fft for each weighted MR signal 
K = zeros(size(img));

for i = 1:size(img,4) 
    thisImg = img(:,:,:,i);
    K(:,:,:,i) = 1/sqrt(length(thisImg(:)))*fftshift(fftn(thisImg));
end