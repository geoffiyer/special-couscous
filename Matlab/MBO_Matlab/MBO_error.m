addpath ImageTools/

if(size(im,3) > 3)
    lidar = im(:,:,4);
    lidar = rescaleIm(lidar);
    im = im(:,:,1:3);
end

im = rescaleIm(im);
K = reshape(dlmread('~/Schoolwork/Chanussot/MBO_Code/Output/classification_result.txt'), [size(im,2) size(im,1)])';
assert(exist('im','var') && exist('lidar','var'),'need image and lidar to check error');

scalings = dlmread('~/Schoolwork/Chanussot/MBO_Code/Output/scaling.txt');

e = segmentationError(cat(3,im/scalings(1),lidar/scalings(2)),K);

rmpath ImageTools/