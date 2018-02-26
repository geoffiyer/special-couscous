function[K] = viewMBO(im)

addpath ../ImageTools/

imSize = [size(im,1) size(im,2) 4];
im = double(permute(reshape(dlmread('~/Schoolwork/Chanussot/MBO_Code/Output/image.txt'),[imSize(2),imSize(1),imSize(3)]),[2 1 3]));
lidar = rescaleIm(im(:,:,4));
im = rescaleIm(im(:,:,1:3));
fidelity = reshape(dlmread('~/Schoolwork/Chanussot/MBO_Code/Output/fidelity.txt'),[size(im,2) size(im,1)])';
fidelity = rescaleIm(fidelity);

K = reshape(dlmread('~/Schoolwork/Chanussot/MBO_Code/Output/classification_result.txt'), [size(im,2) size(im,1)])';

%     try
%         cartoonIm = cartoon(im,K);
%     catch
if(exist('fidelity','var'))
    fidelity = rescaleIm(fidelity);
    lidar = rescaleIm(lidar);
    subplot(2,2,1);
    imshow(im);
    title('RGB image');
    subplot(2,2,2);
    imshow(lidar*64,jet);
    title('lidar image');
    subplot(2,2,3);
    imshow(fidelity*64,jet);
    title('fidelity');
    subplot(2,2,4);
    imshow(rescaleIm(K)*64,jet);
    title('classification');
elseif(exist('lidar','var'))
    subplot(1,3,1);
    imshow(im);
    title('RGB image');
    subplot(1,3,2);
    imshow(lidar);
    title('lidar image');
    subplot(1,3,3);
    imagesc(K);
    title('classification');
else
    subplot(1,2,1);
    imshow(im);
    title('RGB image');
    subplot(1,2,2);
    imagesc(K);
    title('Classification');
end

end