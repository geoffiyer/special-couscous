%% We got some DFC2018 data babay
addpath ../ImageTools/
addpath ../Nystrom/
addpath ../GraphMatch/envi/
addpath ../MBO_Matlab/

%% The variables
percentGT = 1.0;
patch = [101 901 1 601];

%% Names of all the files in question
filedir = '/home/gsiyer/Schoolwork/Chanussot/Matlab/DFC2018_Data/2018_Release_Phase1/';
HSname = 'HSI/2018_IEEE_GRSS_DFC_HSI_TR';
HShdrname = 'HSI/Copy_of_2018_IEEE_GRSS_DFC_HSI_TR.HDR';
GTname = 'GT/2018_IEEE_GRSS_DFC_GT_TR.tif';
lidarnames = [string(strcat('Lidar_GeoTiff_Rasters/DEM_C123_3msr/','UH17_GEG051_TR.tif')) ; ...
              string(strcat('Lidar_GeoTiff_Rasters/DEM_C123_TLI/' ,'UH17_GEG05_TR.tif'))  ; ...
              string(strcat('Lidar_GeoTiff_Rasters/DEM+B_C123/'   ,'UH17_GEM051_TR.tif')) ; ...
              string(strcat('Lidar_GeoTiff_Rasters/DSM_C12/'      ,'UH17c_GEF051_TR.tif')); ...
              string(strcat('Lidar_GeoTiff_Rasters/Intensity_C1/' ,'UH17_GI1F051_TR.tif')); ...
              string(strcat('Lidar_GeoTiff_Rasters/Intensity_C2/' ,'UH17_GI2F051_TR.tif')); ...
              string(strcat('Lidar_GeoTiff_Rasters/Intensity_C3/' ,'UH17_GI3F051_TR.tif'))];
lidarscales = [1; 1; 9e36; 9e36; 2e34; 2e34; 2e34];
lidarnames = cellstr(lidarnames);

%% Load the hyperspectral
[HS,info] = enviread(strcat(filedir,HSname),strcat(filedir,HShdrname));
HS = rescaleIm(double(HS));

rmpath ../GraphMatch/envi/

%% GT means ground truth!!!!
GT = double(imread(strcat(filedir,GTname)));

%% Lots of lidars to investigate
numLidars = 7;

lidars = zeros(size(HS,1),size(HS,2),numLidars);
for i = 1:min(numLidars,size(lidarnames,1))
    lidarorig = rescaleIm(abs(double(imread(strcat(filedir,lidarnames{i})))));
%     lidar = killOutliers(lidarorig);
    lidarorig = rescaleIm(double(lidarorig));
    lidarorig = lidarorig(1:2:size(lidarorig,1),1:2:size(lidarorig,2));
    lidarorig = lidarorig * lidarscales(i);
    lidarorig(lidarorig > 1) = 1;
    lidars(:,:,i) = lidarorig;
%     lidar = rescaleIm(double(lidar));
%     imwrite(lidarorig*64,jet,strcat(writeToFolder,'lidarorig',int2str(i),'.png'));
%     imwrite(lidar*64,jet,strcat(writeToFolder,'lidar',int2str(i),'.png'));
end

% for i = 1:numLidars
%     temp = rescaleIm(double(imread(strcat(writeToFolder,'lidaredit',int2str(i),'.png'))));
%     lidars(:,:,i) = temp(1:2:size(temp,1),1:2:size(temp,2));
% end
data = cat(3,HS,lidars);

xmin = max(patch(1),1);
xmax = min(patch(2), size(data,2));
ymin = max(patch(3),1);
ymax = min(patch(4), size(data,1));
data = data(ymin:ymax, xmin:xmax, :);
GT = GT(1:2:size(GT,1),1:2:size(GT,2));
GT = GT(ymin:ymax, xmin:xmax);
GT = rescaleClasses(GT,100);
HSim = data(:,:,[48 32 16]);

GT(rand(size(GT)) > percentGT) = 0;

[E,dex] = getMultimodalWeights(data,[size(HS,3) ones(1,numLidars)],2,101);
[V,D] = getEvecs(E,dex,101);

%% Try to MBO it LET'S GOOOOO
numClasses = max(GT(:));
fidelity = rescaleIm(GT);
[K, info] = MBO(V, D, numClasses,fidelity);

K = K+1;
K = matchClasses(K,GT);

subplot(1,2,1)
imshow(K/numClasses*64,jet)
subplot(1,2,2)
imshow(rescaleIm(GT)*64,jet)

% %% Kmeans it already
% rng('shuffle');
% K = eff_kmeans(V(1:20),numClasses,200);
% K = reshape(K,[size(data,1) size(data,2)]);
% K = smoothPixels(K,1);

%% Write stuff to a file
if( false )
    writeToFolder = '/home/gsiyer/Schoolwork/Chanussot/MeetingNotes/2018-02-09/DFC2018/';
    imwrite(HSim*0.8+0.2,strcat(writeToFolder,'Hyperspectral_RGB_Bands.png'));    
    imwrite(rescaleIm(GT)*64,jet,strcat(writeToFolder,'Ground_Truth.png'));
    for i = 1:9
        imwrite(rescaleIm(reshape(V(:,i),[size(data,1) size(data,2)]))*64,jet,strcat(writeToFolder,'evec',int2str(i),'.png'))
    end
    for i = 1:numLidars
        imwrite(rescaleIm(data(:,:,size(HS,3)+i))*64,jet,strcat(writeToFolder,'lidar',int2str(i),'.png'));
    end
    imwrite(K/numClasses*64,jet,strcat(writeToFolder,'MBO_Classification.png'));
end