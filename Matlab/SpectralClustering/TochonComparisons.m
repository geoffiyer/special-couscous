addpath ../ImageTools

A = load('Umbrella/Data.mat');
B = load('TochonData/RGB_depth_results.mat');
seg_braid = B.optimalseg_braid;
temp = unique(seg_braid);
for i = 1:size(temp,1)
    seg_braid(seg_braid == temp(i)) = i;
end
seg_BPT = B.optimalseg_consensusBPT;
temp = unique(seg_braid);
for i = 1:size(temp,1)
    seg_braid(seg_braid == temp(i)) = i;
end
nclasses = 6;
string = 'umbrella';

% A = load('DFC/DFC2015.mat'); % For DFC, nclasses = 5 looks pretty good.
%                              % Should probably think more about how to
%                              % choose number.
% B = load('TochonData/ortho_LiDAR_results.mat');
% nclasses = 8;
% string = 'DFC';

displayGraph = false;
pixelsmoothing = 1;
p = inf;

imwrite(A.im,'~/Schoolwork/Chanussot/Papers/Image_Segmentation/Images/Temp/optical.png');
A.lidar = (A.lidar - min(min(A.lidar))) ./ (max(max(A.lidar)) - min(min(A.lidar))) * 64;
imwrite(A.lidar,jet,'~/Schoolwork/Chanussot/Papers/Image_Segmentation/Images/Temp/lidarColor.png');

[~,I,K,scaling1,scaling2] = SpectralClassification(A.im,A.lidar,nclasses,displayGraph,pixelsmoothing,p);
errorMine = segmentationError(cat(3,A.im/scaling1,A.lidar/scaling2),K,'region');
ccMine = connComp(K);
numRegionsMine = max(max(ccMine));
cartoonMine = cartoon(A.im,K);
imwrite(cartoonMine,strcat('~/Schoolwork/Chanussot/Papers/Image_Segmentation/Images/Temp/', ...
    string, 'CartoonGraph.png'));

errorBraid = segmentationError(cat(3,A.im/scaling1,A.lidar/scaling2),seg_braid,'region');
numRegionsBraid = max(max(seg_braid));
cartoonBraid = cartoon(A.im,seg_braid);
imwrite(cartoonBraid,strcat('~/Schoolwork/Chanussot/Papers/Image_Segmentation/Images/Temp/', ...
    string, 'CartoonBraid.png'));

errorBPT = segmentationError(cat(3,A.im/scaling1,A.lidar/scaling2),seg_BPT,'region');
numRegionsBPT = max(max(seg_BPT));
cartoonBPT = cartoon(A.im,seg_BPT);
imwrite(cartoonBPT,strcat('~/Schoolwork/Chanussot/Papers/Image_Segmentation/Images/Temp/', ...
    string, 'CartoonBPT.png'));

subplot(1,3,1)
imshow(rescaleIm(K)*64,jet)
subplot(1,3,2)
imshow(rescaleIm(seg_braid)*64,jet)
subplot(1,3,3)
imshow(rescaleIm(seg_BPT)*64,jet)

% imwrite(rescaleIm(K)*64,jet,'/home/gsiyer/Schoolwork/Chanussot/MeetingNotes/2018-01-18/specClust.png')
% imwrite(rescaleIm(seg_braid)*64,jet,'/home/gsiyer/Schoolwork/Chanussot/MeetingNotes/2018-01-18/braid.png')
% imwrite(rescaleIm(seg_BPT)*64,jet,'/home/gsiyer/Schoolwork/Chanussot/MeetingNotes/2018-01-18/BPT.png')