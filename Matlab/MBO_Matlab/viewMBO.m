function[K] = viewMBO(im)

addpath ../ImageTools/

imSize = [size(im,1) size(im,2) size(im,3)];
lidar = im(:,:,4);
im = im(:,:,1:3);
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
    imshow((K+1)/(size(unique(fidelity),1))*64,jet);
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