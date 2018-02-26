    
%% Load the data

addpath ../ImageTools

useNLM = false;
useHungarian = false;
normType = 'L2';
removeEdgeweightOutliers = false;
clusterSize = 100;
nonlinear = true;
superSize = 1000;

% Edurne's Change detection beach dunes pictures: December 2017
% Something was done to the registration (I didn't ask any questions)
imgpath = '/home/gsiyer/Schoolwork/Chanussot/Edurne/MaspalomasSubsets/MASP_2013aug11.dat';
hdrpath = '/home/gsiyer/Schoolwork/Chanussot/Edurne/MaspalomasSubsets/MASP_2013aug11.hdr';
addpath('envi');
[pic1,~] = enviread(imgpath,hdrpath);
rmpath('envi');
imgpath = '/home/gsiyer/Schoolwork/Chanussot/Edurne/MaspalomasSubsets/MASP_2015jun04.dat';
hdrpath = '/home/gsiyer/Schoolwork/Chanussot/Edurne/MaspalomasSubsets/MASP_2015jun04.hdr';
addpath('envi');
[pic2,~] = enviread(imgpath,hdrpath);
rmpath('envi');
clear imgpath hdrpath
pic1 = rescaleIm(pic1);
pic2 = rescaleIm(pic2);

pic1 = pic1(250:620,600:1101,:);
pic2 = pic2(250:620,600:1101,:);

% % Edurne's Change detection beach dunes pictures: November 2017
% pic1 = load('/home/gsiyer/Schoolwork/Chanussot/Edurne/MaspalomasSubsets/MASP_2013aug11.mat');
% pic1 = double(pic1.I);
% pic1 = pic1(1:500,:,:);
% pic2 = load('/home/gsiyer/Schoolwork/Chanussot/Edurne/MaspalomasSubsets/MASP_2015jun04.mat');
% pic2 = double(pic2.I);
% pic2 = pic2(1:500,:,:);
% pic1 = pic1(:,:,[5 3 2]);
% pic2 = pic2(:,:,[5 3 2]);
% pic1 = pic1/max(max(max(pic1)));
% pic2 = pic2/max(max(max(pic2)));

% % Test the nonlinear alteration plan
% patch = [801 1000 1601 1800];
% pic1 = load_data('SPOT/','W_Spot-Before',patch);
% pic2 = load_data('SPOT/','W_Spot-Before',patch);
% nonlinear = true;

% % New idea: load both from before-flood (not 1 before 1 after), but a bit offset
% patch = [801 1000 1601 1800];
% offset = 20;
% pic1 = load_data('SPOT/','W_Spot-Before',patch);
% pic2 = load_data('SPOT/','W_Spot-Before',patch + [offset offset offset offset]);

% % Try loading in the RGB-Lidar stuff from the image segmentation paper.
% temp = load('../SpectralClustering/DFC/DFC2015.mat');
% pic1 = temp.im;
% pic2 = temp.lidar;
% clear temp;

% % Try loading in the RGB-Lidar stuff from the image segmentation paper.
% temp = load('../SpectralClustering/Umbrella/Data.mat');
% pic1 = rescaleIm(temp.im);
% pic2 = rescaleIm(temp.lidar);
% clear temp;

% % Try loading in the RGB-Lidar stuff from the image segmentation paper.
% temp = load('../SpectralClustering/Jadeplant/Data.mat');
% pic1 = rescaleIm(temp.im);
% pic2 = rescaleIm(temp.lidar);
% clear temp;
 
% % The hyperspectral stuff
% temp = load('../../HyperspectralData/SWIR.mat');
% patch = [1 271 1 271];
% pic1 = rescaleIm(temp.M(patch(1):patch(2),patch(3):patch(4),4:6));
% pic2 = rescaleIm(temp.M(patch(1):patch(2),patch(3):patch(4),124:126));
% clear temp patch;

% % Flood data: algorithm fails. Somehow eigenvector signs are mixed up?
% patchx = 1201;
% patchy = 1401;
% patch = [patchx patchx+299 patchy patchy+299];
% pic1 = load_data('SPOT/','W_Spot-Before',patch);
% pic2 = load_data('SPOT/','W_Spot-After',patch);

% % Flood data: similar before/after in this patch. 2 choices
% patch = [1201 1700 1701 2100];
% % patch = [1001 1250 1401 1650];
% pic1 = load_data('SPOT/','W_Spot-Before',patch);
% temp1 = max(max(pic1));
% pic2 = load_data('SPOT/','W_Spot-After',patch);
% temp2 = max(max(pic2));
% temp3 = max(temp1,temp2);
% pic1 = pic1./temp3;
% pic2 = pic2./temp3;
% % pic1 = load_data('ERS/','W_23538-Before',patch);
% % pic2 = load_data('ERS/','W_29049-After',patch);
% clear temp1 temp2 temp3;

% % Flood data: different before/after in this patch
% patch = [851 950 2551 2650];
% pic1 = load_data('SPOT/','W_Spot-Before',patch);
% pic2 = load_data('SPOT/','W_Spot-After',patch);

% % Synthetic Data: The spot thingy
% [pic1,pic2] = GenerateData(6,200,200);
% super = sub2ind([200 200],[1 200 200 1],[1 1 200 200]);

% % Synthetic Data: Normal picture with lots of small changes and 1 big.
% temp = load('../SpectralClustering/Adirondack/Data.mat');
% pic1 = temp.im;
% pic2 = pic1;
% pic2(20:45,20:45,:) = 0;
% clear temp;

assert(isequal(size(pic1(1:2)),size(pic2(1:2))),'pic1 and pic2 must be corregistered');

% try adding nlm
if(useNLM)
    nlmpic1 = nonlocal(pic1,1,1);
    nlmpic2 = nonlocal(pic2,1,1);
else
    nlmpic1 = pic1;
    nlmpic2 = pic2;
end

pic1Size = [size(nlmpic1,1),size(nlmpic1,2),size(nlmpic1,3)];
pic2Size = [size(nlmpic2,1),size(nlmpic2,2),size(nlmpic2,3)];
temp = randi(pic1Size(1)*pic1Size(2),superSize,1);
super = [temp, temp];

%% This was editing super in one VERY specific case.
% i = 1;
% while( i <= size(super,1))
%     [temp1 temp2] = ind2sub([pic1Size(1), pic1Size(2)], super(i,1));
%     if( 20 <= temp1 && temp1 <= 50)
%         super(i,:) = [];
%         i = i-1;
%     elseif( 20 <= temp2 && temp2 <= 50)
%         super(i,:) = [];
%         i = i-1;
%     end
%     i = i+1;
% end
% clear i temp1 temp2;

%% Change the data into a n x 3 vector
X1 = reshape(nlmpic1, [pic1Size(1)*pic1Size(2) pic1Size(3)]);
X2 = reshape(nlmpic2, [pic2Size(1)*pic2Size(2) pic2Size(3)]);
if(nonlinear)
    X2 = nonlinearAlteration(X2, pic2Size);
end

%% Show what we loaded
if( pic1Size(3) <= 3)
    subplot(2,3,1)
    imshow(reshape(X1,pic1Size));
    title('Picture X');
end
if( pic2Size(3) <= 3)
    subplot(2,3,4)
    imshow(reshape(X2,pic2Size));
    title('Picture Y');
end

[U1,U2] = getEigenvectors(X1,X2,normType);
[signs, evecassign] = fixSigns(U1(super(:,1),:),U2(super(:,2),:));
U2 = U2.*signs';
U2(:, evecassign(:,2)) = U2;

%% Do the graph match
if(useHungarian)
    [assign,Y1,Y2,Yassign,U1,U2] = GraphMatch(X1,X2,super,normType,clusterSize);
    assign1to2 = sortrows(assign,1);
    assign2to1 = sortrows(assign,2);
    clear assign;
else
    [assign1to2,assign2to1] = GraphMatch2(U1,U2,clusterSize);
end

% weights corresponding to the edges we chose, as well as a threshold.
edgeweights1to2 = dot(U1(assign1to2(:,1),:),U2(assign1to2(:,2),:),2);
% edgethresh1to2 = FindEdgethresh(edgeweights1to2);
edgeweights2to1 = dot(U1(assign2to1(:,1),:),U2(assign2to1(:,2),:),2);
% edgethresh2to1 = FindEdgethresh(edgeweights2to1);

% %% Look at edge weights a little
% % this code picks out the strongNum of the strongest edges, and gives you the assignment
% strongNum = 30;
% edgeweights = [assign(:,1), assign(:,2), U(sub2ind(size(U), assign(:,2), assign(:,1)))];
% edgeweights = sortrows(edgeweights,-3);
% strongAssign = edgeweights(1:strongNum,1:2);

%% Difference between coregister and graph match
X1Norms = sqrt(sum(X1.^2,2));
X2Norms = sqrt(sum(X2.^2,2));
diffNorms = zeros(size(assign1to2,1),4);

% Calculate difference as ratio of (perm1 - perm2) / (perm1 + perm2)
% Idea: We need some sort of normalization
diff_vec = sqrt(2)*(X1(assign1to2(:,2),:) - X1(assign1to2(:,1),:)) ./ ...
    (X1Norms(assign1to2(:,2)) + X1Norms(assign1to2(:,1)));
diffNorms(:,1) = sqrt(sum(diff_vec.^2,2));

diff_vec = sqrt(2)*(X2(assign1to2(:,2),:) - X2(assign1to2(:,1),:)) ./ ...
    (X2Norms(assign1to2(:,2)) + X2Norms(assign1to2(:,1)));
diffNorms(:,2) = sqrt(sum(diff_vec.^2,2));

% % We're focusing on pairs that GraphMatch says are good, but diffNorms says
% % are bad. This sort of makes sense but I'm still working on understanding.
% diffNorms(edgeweights1to2<edgethresh1to2,1:2) = ... 
%     repmat([0 0], [size(edgeweights1to2(edgeweights1to2 < edgethresh1to2)) , 1]);

diff_vec = sqrt(2)*(X1(assign2to1(:,2),:) - X1(assign2to1(:,1),:)) ./ ...
    (X1Norms(assign2to1(:,2)) + X1Norms(assign2to1(:,1)));
diffNorms(:,3) = sqrt(sum(diff_vec.^2,2));

diff_vec = sqrt(2)*(X2(assign2to1(:,2),:) - X2(assign2to1(:,1),:)) ./ ...
    (X2Norms(assign2to1(:,2)) + X2Norms(assign2to1(:,1)));
diffNorms(:,4) = sqrt(sum(diff_vec.^2,2));

% diffNorms(edgeweights2to1<edgethresh2to1,3:4) = ...
%     repmat([0 0], [size(edgeweights2to1(edgeweights2to1 < edgethresh2to1)) , 1]);

%% alternate version of diffnorms where you look at the exact difference rather than some ratio
% diff_vec = X1(assign(:,2),:) - X1(assign(:,1),:);
% diffNorms(:,1) = sqrt(sum(diff_vec.^2,2));
% 
% diff_vec = X2(assign(:,2),:) - X2(assign(:,1),:);
% diffNorms(:,2) = sqrt(sum(diff_vec.^2,2));
% 
% assign = sortrows(assign,2);
% diff_vec = X1(assign(:,2),:) - X1(assign(:,1),:);
% diffNorms(:,3) = sqrt(sum(diff_vec.^2,2));
% 
% diff_vec = X2(assign(:,2),:) - X2(assign(:,1),:);
% diffNorms(:,4) = sqrt(sum(diff_vec.^2,2));

%% Edgeweight outliers

if(removeEdgeweightOutliers)
addpath ../Nystrom
[~,centers,~] = eff_kmeans(edgeweights1to2,2,100);
centers = sort(centers);
edgeweights1to2(edgeweights1to2 < centers(1)) = centers(1);
[~,centers,~] = eff_kmeans(edgeweights2to1,2,100);
centers = sort(centers);
edgeweights2to1(edgeweights2to1 < centers(1)) = centers(1);
rmpath ../Nystrom/
clear centers
end

clear removeEdgeweightOutliers
%% Show the edited

if( pic2Size(3) <= 3)
    subplot(2,3,2)
    imshow(reshape(X2(assign1to2(:,2),:), pic2Size));
    title('perm(Y)');
end

if( pic1Size(3) <= 3)
    subplot(2,3,5)
    imshow(reshape(X1(assign2to1(:,1),:), pic1Size));
    title('perm^{-1}(X)');
end

% subplot(2,4,3)
% imshow(64*rescaleIm(reshape(edgeweights1to2,[pic1Size(1),pic1Size(2)])),jet)
% title('edgeweights, the other')

diffIm1 = rescaleIm(reshape(diffNorms(:,2), [pic1Size(1),pic1Size(2)]));
subplot(2,3,3)
imshow(64*diffIm1, jet);
title('norm( Y - perm(Y) )');

% subplot(2,4,7)
% imshow(64*rescaleIm(reshape(edgeweights2to1,[pic1Size(1),pic1Size(2)])),jet)
% title('edgeweights, one permutation')

diffIm2 = rescaleIm(reshape(diffNorms(:,3), [pic1Size(1),pic1Size(2)]));
subplot(2,3,6)
imshow(64*diffIm2, jet);
title('norm( X - perm^{-1}(X) )');

%% This is all plotting specific things for images for slides and stuff

% naiveDiff = reshape( sqrt(sum( (X1-X2).^2,2)), [size(pic1,1) size(pic1,2)]);
% imshow(naiveDiff)
% imwrite(naiveDiff,'~/Schoolwork/Chanussot/ATC/Images/GraphMatch/naivediff.png');
% imwrite(64*reshape(diffNorms(:,3), [size(pic1,1), size(pic1,2)]), jet, ...
%     '~/Schoolwork/Chanussot/ATC/Images/GraphMatch/normMap.png');
% imwrite(pic1,'~/Schoolwork/Chanussot/ATC/Images/GraphMatch/dataBefore.png');
% imwrite(pic2,'~/Schoolwork/Chanussot/ATC/Images/GraphMatch/dataAfter.png');

% figure
% assign = sortrows(assign,1);
% perm = assign(:,2);
% invperm(perm) = 1:size(assign,1);
% subplot(2,3,1)
% imshow(reshape(X1,size(pic1)))
% subplot(2,3,2)
% imshow(reshape(X1(perm,:),size(pic1)))
% subplot(2,3,3)
% imshow(reshape(X1(invperm,:),size(pic1)))
% subplot(2,3,4)
% imshow(reshape(X2,size(pic1)))
% subplot(2,3,5)
% imshow(reshape(X2(perm,:),size(pic1)))
% subplot(2,3,6)
% imshow(reshape(X2(invperm,:),size(pic1)))

% subplot(1,3,1)
% imshow(pic2)
% title('Picture X (before change)')
% subplot(1,3,2)
% imshow(pic1)
% title('Picture Y (after change)')
% subplot(1,3,3)
% imshow(64*reshape(diffNorms(:,4), [size(pic1,1), size(pic1,2)]), jet);
% title('norm( X(i) - X(rho(i) )');

% imwrite(pic1,'~/Schoolwork/Chanussot/Papers/Research_Statement_06-2017/Images/ChangeDetection/Example2/dataBefore.png','png');
% imwrite(pic2,'~/Schoolwork/Chanussot/Papers/Research_Statement_06-2017/Images/ChangeDetection/Example2/dataAfter.png','png');
% imwrite(64*reshape(diffNorms(:,4), [size(pic1,1), size(pic1,2)]), jet, ...
%     '~/Schoolwork/Chanussot/Papers/Research_Statement_06-2017/Images/ChangeDetection/Example2/normMap.png','png');

% figure
% subplot(2,2,1)
% imshow(pic1)
% title('Before Flood')
% subplot(2,2,3)
% imshow(pic2)
% title('After Flood')
% subplot(2,2,4)
% imshow(64*rescaleIm(reshape(diffNorms(:,2), [pic1Size(1),pic1Size(2)])), jet);
% title('Norm of (perm_{X to Y}(Y) - Y)')
% subplot(2,2,2)
% imshow(64*rescaleIm(reshape(diffNorms(:,3), [pic1Size(1),pic1Size(2)])), jet);
% title('Norm of (perm_{Y to X}(X) - X)')

clear useNLM X1Norms X2Norms diff_vec patch picSize perm invperm offset ...
    patchx patchy tempedgeweights useHungarian normType clusterSize