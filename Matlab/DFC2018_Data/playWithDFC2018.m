%% We got some DFC2018 data babay
addpath ../ImageTools/
addpath ../Nystrom/
addpath ../GraphMatch/envi/
addpath ../MBO_Matlab/

%% The variables
percentGT = 1.0;
patch = [1 300 1 601];
norms = [0 0 0 0 0 0 0 0];
pNorm = 2;
nEvecs = 201;
fidelityConst = 1.0e3;
dtConst = 1.0e-1;
calculateEvecs = true;
loadfiles = true;

if(loadfiles)
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
    [HS,~] = enviread(strcat(filedir,HSname),strcat(filedir,HShdrname));
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
end

xmin = max(patch(1),1);
xmax = min(patch(2), size(data,2));
ymin = max(patch(3),1);
ymax = min(patch(4), size(data,1));
data = data(ymin:ymax, xmin:xmax, :);
HSim = data(:,:,[48 32 16]);
totalNumClasses = max(GT(:));
GT = GT(1:2:size(GT,1),1:2:size(GT,2));
GT = GT(ymin:ymax, xmin:xmax);

% This comes from looking at the actual data.
% Some classes are indistinguishable or at least borderline
% GT(GT == 2) = 1; % 2 and 1 are both roofs of buldings. Not sure exactly which.
% GT(GT == 5) = 4; % 4 and 5 are both grass types I think??
% GT(GT == 3) = 13; % 3, 10, 11, 12, 13, 16 are all types of road
% GT(GT == 10) = 13; 
% GT(GT == 11) = 13; 
% GT(GT == 12) = 13;
% GT(GT == 14) = 13;
% GT(GT == 16) = 13;

[GTedit,classesIdx] = rescaleFidelity(GT); % if GTedit(i,j) == k, then GT(i,j) == classesIdx(k)
numClasses = max(GTedit(:));
GTedit(rand(size(GTedit)) > percentGT) = 0;

%% pick index for nystrom based on the semisupervised classes
dex = [];
% classSizes = zeros(numClasses,1);
% for i = 1:numClasses
%     classSizes(i) = sum(GT == i);
% end
amt = ceil(nEvecs/numClasses);
for i = 1:numClasses
    temp = find(GTedit == i);
    perm = randperm(size(temp,1));
    dex = [dex;temp(perm( 1:min(amt,size(temp,1)) ))];
end
if(size(dex,1) < nEvecs)
    temp = find(GTedit == 0);
    perm = randperm(size(temp,1));
    amt = nEvecs - size(dex,1);
    dex = [dex;temp(perm( 1:min(amt,size(temp,1)) ))];
end
dex = dex(randperm(size(dex,1)));

if(calculateEvecs)
    [E,dex] = getMultimodalWeights(data,[size(HS,3) ones(1,numLidars)],pNorm,nEvecs,norms,dex);
    [V,D] = getEvecs(E,dex,nEvecs);
    V = orth(V);
end

%% Try to MBO it LET'S GOOOOO
fidelity = GTedit/255;
[K_old, info] = MBO(V, D, numClasses,fidelity,fidelityConst,dtConst);

K_old = K_old+1;
K = zeros(size(K_old));

% change class labeling from MBO to original GT
for i = 1:size(classesIdx,1)
    K(K_old == i) = classesIdx(i);
end

confus = confusionMat(K,GT); % confus(i,j) = #points classified as i, and ground truth as j
confus = confus./sum(confus,1);

data = reshape(data, [size(data,1)*size(data,2) size(data,3)]);
GTvec = reshape(GT, [size(GT,1)*size(GT,2) size(GT,3)]);

dataMeans = zeros(numClasses,size(data,2));

for i = 1:max(GT(:))
    dataMeans(i,:) = mean(data(GTvec == i,:));
end

% Differences between classes.
% Use this to get an idea of where to combine
A = 1-real(angledist(dataMeans,dataMeans));

subplot(1,2,1)
imshow(K/totalNumClasses*64,jet)
subplot(1,2,2)
imshow(GT/totalNumClasses*64,jet)

% %% Kmeans it already
% rng('shuffle');
% K = eff_kmeans(V(1:20),numClasses,200);
% K = reshape(K,[size(data,1) size(data,2)]);
% K = smoothPixels(K,1);

%% Write stuff to a file
if( false )
    writeToFolder = '/home/gsiyer/Schoolwork/Chanussot/MeetingNotes/2018-03-12/DFC2018/';
    imwrite(HSim*0.8+0.2,strcat(writeToFolder,'Hyperspectral_RGB_Bands.png'));    
    imwrite(GT/totalNumClasses*64,jet,strcat(writeToFolder,'Ground_Truth.png'));
    data = reshape(data,[size(HSim,1) size(HSim,2) 57]);
    for i = 1:9
        imwrite(rescaleIm(reshape(V(:,i),[size(HSim,1) size(HSim,2)]))*64,jet,strcat(writeToFolder,'evec',int2str(i),'.png'))
    end
    for i = 1:numLidars
        imwrite(rescaleIm(data(:,:,size(HS,3)+i))*64,jet,strcat(writeToFolder,'lidar',int2str(i),'.png'));
    end
    imwrite(K/totalNumClasses*64,jet,strcat(writeToFolder,'MBO_Classification.png'));
    for i = 1:size(dataMeans,1)
        if(~isnan(dataMeans(i,1)))
            plot(dataMeans(i,:));
            print(strcat(writeToFolder,'Class_',int2str(i),'_mean.png'),'-dpng');
        end
    end
    confus = [confus, zeros(size(confus,1),totalNumClasses-size(confus,2))];
    confus = [confus; zeros(totalNumClasses - size(confus,1),size(confus,2))];
    imagesc(confus)
    print(strcat(writeToFolder,'confusionMatrix.png'),'-dpng');
    A = [A, zeros(size(A,1),totalNumClasses-size(A,2))];
    A = [A; zeros(totalNumClasses - size(A,1),size(A,2))];
    temp = min(A(A ~= 0));
    A(A == 0) = temp;
    imagesc(A)
    print(strcat(writeToFolder,'spectralComparisons.png'),'-dpng');
end