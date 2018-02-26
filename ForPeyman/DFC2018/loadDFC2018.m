%% Load DFC2018 data
addpath ../envi/

%% Names of all the files in question

% Peyman: You will have to change this to the directory where you have the
%     DFC2018 data.
filedir = '/home/gsiyer/Schoolwork/Chanussot/Matlab/DFC2018_Data/2018_Release_Phase1/';

% Names of individual files in dataset. I expect this doesn't need to be
% changed.
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
lidarnames = cellstr(lidarnames);

%% Load the hyperspectral into HS
[HS,info] = enviread(strcat(filedir,HSname),strcat(filedir,HShdrname));

rmpath ../envi/

%% Load the ground truth into GT
GT = imread(strcat(filedir,GTname));
GT = rescaleIm(double(GT));

%% 7 different lidar images to load
for i = 1:size(lidarnames,1)
    lidarorig = rescaleIm(abs(double(imread(strcat(filedir,lidarnames{i})))));
    lidar = killOutliers(lidarorig);
    lidarorig = rescaleIm(double(lidarorig));
    lidar = rescaleIm(double(lidar));
end


numLidars = 7;
lidars = zeros(size(HS,1),size(HS,2),numLidars);
for i = 1:numLidars
    temp = rescaleIm(double(imread(strcat(writeToFolder,'lidaredit',int2str(i),'.png'))));
    lidars(:,:,i) = temp(1:2:size(temp,1),1:2:size(temp,2));
end
data = cat(3,HS,lidars);

xmin = max(patch(1),1);
xmax = min(patch(2), size(data,2));
ymin = max(patch(3),1);
ymax = min(patch(4), size(data,1));
data = data(ymin:ymax, xmin:xmax, :);
GT = GT(1:2:size(GT,1),1:2:size(GT,2));
GT = GT(ymin:ymax, xmin:xmax);
HSim = data(:,:,[48 32 16]);

[E,dex] = getMultimodalWeights(data,[size(HS,3) ones(1,numLidars)],10,101);
[V,D] = getEvecs(E,dex,101);

%% Kmeans it already
rng('shuffle');
K = eff_kmeans(V,numClasses,200);
K = reshape(K,[size(data,1) size(data,2)]);
K = smoothPixels(K,1);

%% Write stuff to a file
if( false )
    imwrite(HSim*0.8+0.2,strcat(writeToFolder,'Hyperspectral.png'));    
    imwrite(GT*64,jet,strcat(writeToFolder,'Ground_Truth.png'));
    for i = 1:9
        imwrite(rescaleIm(reshape(V(:,i),[601 601]))*64,jet,strcat(writeToFolder,'evec',int2str(i),'.png'))
    end
    for i = 1:numLidars
        imwrite(rescaleIm(data(:,:,size(HS,3)+i))*64,jet,strcat(writeToFolder,'lidar',int2str(i),'.png'));
    end
    imwrite(rescaleIm(K)*64,jet,strcat(writeToFolder,'specClust.png'));
end