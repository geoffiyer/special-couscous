
addpath ../ImageTools/

displayCartoons = false;

% A = load('DFC/DFC2013.mat');
% baseIm = A.im(:,:,[40,70,100]);
% unionClasses1 = 3;
% unionClasses2 = 3;
% nclasses = 6;

% A = load('Jadeplant/Data.mat');
% unionClasses1 = 3;
% unionClasses2 = 2;
% nclasses = 7;

A = load('Umbrella/Data.mat');
unionClasses1 = 3;
unionClasses2 = 2;
nclasses = 6;

% A = load('DFC/DFC2015.mat'); % For DFC, nclasses = 5 looks pretty good.
%                              % Should probably think more about how to
%                              % choose number.
% unionClasses1 = 3;
% unionClasses2 = 2;
% nclasses = 6;

displayGraph = false;
pixelsmoothing = 1;
p = inf;

imwrite(A.im,'~/Schoolwork/Chanussot/Papers/Image_Segmentation/Images/Temp/optical.png');
A.lidar = (A.lidar - min(min(A.lidar))) ./ (max(max(A.lidar)) - min(min(A.lidar))) * 64;
imwrite(A.lidar,jet,'~/Schoolwork/Chanussot/Papers/Image_Segmentation/Images/Temp/lidarColor.png');

tic
[errorMine,I,K,scaling1,scaling2] = SpectralClassification(A.im,A.lidar,nclasses,displayGraph,pixelsmoothing,p);
cartoonMine = cartoon(A.im,K);
imwrite(rescaleIm(K)*64,jet,'~/Schoolwork/Chanussot/Papers/Image_Segmentation/Images/Temp/specClust.png');
timeMine = toc;

tic
[~,~,tempK1] = SpectralClassification(A.im,A.im,unionClasses1,displayGraph,pixelsmoothing,p);
[~,~,tempK2] = SpectralClassification(A.lidar,A.lidar,unionClasses2,displayGraph,pixelsmoothing,p);
unionK = unionClasses1*(tempK2-1) + tempK1;
cartoonUnion = cartoon(A.im,unionK);
errorUnion = segmentationError(cat(3,A.im/scaling1,A.lidar/scaling2), unionK, 'regions');
% ...              / (sqrt(size(A.im,3)+size(A.lidar,3)) * size(A.im,1) * size(A.im,2));

imWithBordersUnion = borders(A.im,unionK);
imwrite(rescaleIm(unionK)*64,jet,'~/Schoolwork/Chanussot/Papers/Image_Segmentation/Images/Temp/intersection.png');
timeUnion = toc;

tic
[~,~,K] = SpectralClassification(A.im,A.im,nclasses,displayGraph,pixelsmoothing,p);
errorOptical = segmentationError(cat(3,A.im/scaling1,A.lidar/scaling2), K, 'regions');
cartoonOptical = cartoon(A.im,K);
imwrite(rescaleIm(K)*64,jet,'~/Schoolwork/Chanussot/Papers/Image_Segmentation/Images/Temp/opticalOnly.png');
timeOptical = toc;

tic
[~,~,K] = SpectralClassification(A.lidar,A.lidar,nclasses,displayGraph,pixelsmoothing,p);
errorLidar = segmentationError(cat(3,A.im/scaling1,A.lidar/scaling2), K, 'regions');
cartoonLidar = cartoon(A.im,K);
imWithBordersLidar = borders(A.im,K);
imwrite(imWithBordersLidar,'~/Schoolwork/Chanussot/Papers/Image_Segmentation/Images/Temp/lidarOnly.png');
timeLidar = toc;

tic
[error2Norm,~,K2Norm] = SpectralClassification(A.im,A.lidar,nclasses,displayGraph,pixelsmoothing,2);
cartoon2Norm = cartoon(A.im,K2Norm);
imwrite(imWithBorders2Norm,'~/Schoolwork/Chanussot/Papers/Image_Segmentation/Images/Temp/2norm.png');
time2Norm = toc;

tic
A.im = rescaleIm(A.im)*1.9;
A.lidar = rescaleIm(A.lidar);
data = cat(3,A.im,A.lidar);
dataReshaped = reshape(data,[size(data,1)*size(data,2) size(data,3)]);
idx = eff_kmeans(dataReshaped, nclasses, 200);
K = reshape(idx,[size(data,1) size(data,2)]);
%imagesc(K);
errorKmeans = segmentationError(cat(3,A.im/scaling1,A.lidar/scaling2), K, 'regions');
imWithBordersKmeans = borders(A.im,K);
cartoonKmeans = cartoon(A.im,K);
imwrite(rescaleIm(K)*64,jet,'~/Schoolwork/Chanussot/Papers/Image_Segmentation/Images/Temp/kmeans.png');
timeKmeans = toc;

h = figure;
h.Position = [300,100,1200,850];
if(displayCartoons)
    subplot(2,3,1)
    imshow(cartoonMine);
    title('ours');
    subplot(2,3,2)
    imshow(cartoonUnion);
    title('Product');
    subplot(2,3,3)
    imshow(cartoonKmeans);
    title('kmeans');
    subplot(2,3,4)
    imshow(cartoon2Norm);
    title('2norm');
    subplot(2,3,5)
    imshow(cartoonOptical);
    title('optical');
    subplot(2,3,6)
    imshow(cartoonLidar);
    title('lidar');    
else
    subplot(2,3,1)
    imshow(imWithBorders);
    title('ours');
    subplot(2,3,2)
    imshow(imWithBordersUnion);
    title('Intersection');
    subplot(2,3,3)
    imshow(imWithBordersKmeans);
    title('kmeans');
    subplot(2,3,4)
    imshow(imWithBorders2Norm);
    title('2norm');
    subplot(2,3,5)
    imshow(imWithBordersOptical);
    title('optical');
    subplot(2,3,6)
    imshow(imWithBordersLidar);
    title('lidar');
end