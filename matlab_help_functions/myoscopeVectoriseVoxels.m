function data = myoscopeVectoriseVoxels(img, mask)
% VECTORISE_PIXELS constitues the third step in the biophysical modeling
% pipeline 
%                                           
% input:
%        path2Voxels      :           File path to save voxels
%        imageFolder      :           Folder path to reconstructed MR scans
%        path2mask        :           File path to the binary mask
%        sliceNo          :           The slice number
%        slicePlane       :           'xz', 'yz', or 'xy'. The plane of
%                                     analysis with standard (y, x, z)
%        fullHeartAnalysis:           If ture, sliceNo and slicePlane are
%                                     ignored.
%                      
%                                                                         
% Author: Mohsen Farzi                                                    
% Email : m.farzi@leeds.ac.uk                                             
% \\

imgDim = size(img);
if length(imgDim)==3
    imgDim = [imgDim, 1];
end

if nargin == 1
    mask = true(imgDim(1:3));
end

idx = find(mask);
[iI, iJ, iK] = ind2sub(imgDim(1:3), idx);

img = reshape(img, [prod(imgDim(1:3)), imgDim(4)]);
data = [iI, iJ, iK, img(idx,:)];

