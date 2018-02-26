%% Takes two individual classifications and products them, then compares with multimodal method.
% Yeah I realize that "union.m" is a dumb name for the file. Too late now.

% A = load('Umbrella/Data.mat');
A = load('DFC/DFC2015.mat');
unionClasses1 = 2;
unionClasses2 = 3;


tic
[myError,~,K,~,imWithBorders,scaling1,scaling2] = classify(A.im,A.lidar,unionClasses1*unionClasses2);
time1 = toc;

tic
[~,~,tempK1,bw1] = classify(A.im,A.im,unionClasses1);
[~,~,tempK2,bw2] = classify(A.lidar,A.lidar,unionClasses2);

unionK = unionClasses1*(tempK2-1) + tempK1;

unionError = segmentationError( cat(3,A.im/scaling1,A.lidar/scaling2), unionK);
% ...              / (sqrt(size(A.im,3)+size(A.lidar,3)) * size(A.im,1) * size(A.im,2));

[~,unionImWithBorders] = borders(A.im,unionK);
time2 = toc;

h = figure;
h.Position = [300,100,1200,850];
subplot(1,2,1)
imshow(imWithBorders);
title('Classification with multimodal method');
subplot(1,2,2)
imshow(unionImWithBorders);
title('Union of individual classifications');

%
%
%

% %% I forgot what this is. I think it's some error-finding code.
% 
% displayGraph = false;
% pixelsmoothing = 2;
% p = inf;
% 
% % Umbrella:
% % A = load('Umbrella/Data.mat');
% % nclasses = 8;
% 
% % DFC:
% A = load('DFC/DFC2015.mat');
% nclasses = 6;
% 
% %Run the code on the dataset of choice
% [~,~,K,~,imWithBorders,~,~,~,~] = classify(A.im,A.im,nclasses,displayGraph,pixelsmoothing,p);
% opticalError = segmentationError(cat(3,A.im/scaling1,A.lidar/scaling2), K);
% imwrite(imWithBorders,'~/Schoolwork/Chanussot/Papers/Image_Segmentation/Images/DFC2015/opticalOnly.png');
% 
% [~,I,K,bw,imWithBorders,~,~,E1,E2] = classify(A.lidar,A.lidar,nclasses,displayGraph,pixelsmoothing,p);
% lidarError = segmentationError(cat(3,A.im/scaling1,A.lidar/scaling2), K);
% 
% %
% %
% %
% 
% A = load('Umbrella/Data.mat');
% data = cat(3,A.im,A.lidar);
% dataReshaped = reshape(data,[size(data,1)*size(data,2) size(data,3)]);
% idx = eff_kmeans(dataReshaped, nclasses, 200);
% K = reshape(idx,[size(data,1) size(data,2)]);
% imagesc(K);
% kmeansError = segmentationError(cat(3,A.im/scaling1,A.lidar/scaling2), K);
