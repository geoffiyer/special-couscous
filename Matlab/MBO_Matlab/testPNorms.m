
addpath ../Nystrom
addpath ../ImageTools

%% Initial input variables stuff whoopie
nystromSize = 100;
nEvecs = 50;
numClasses = 6;

% % The data
% A = load('../SpectralClustering/Umbrella/Data.mat');
% im = A.im;
% lidar = A.lidar;
% lidar = rescaleIm(lidar);
% im = rescaleIm(im);

im = zeros(150,150,3);
lidar = zeros(150,150,3);

for i = 1:size(im,1)
    for j = 1:size(im,2)
        if i <= size(im,1)/2
            im(i,j,:) = [1 1 0];
        else
            im(i,j,:) = [0 1 1];
        end
        if j <= size(im,2)/2
            lidar(i,j,:) = [1 1 0];
        else
            lidar(i,j,:) = [0 1 1];
        end
    end
end
im = im + randn(size(im,1),size(im,2))*0.002;
lidar = lidar + randn(size(im,1),size(im,2))*0.002;

fidelity = zeros([size(im,1) size(im,2)]);
fidelity(40:45,40:45) = 1/255;
fidelity(100:105,100:105) = 2/255;
fidelity(40:45,100:105) = 3/255;
fidelity(100:105,40:45) = 4/255;

subplot(1,2,1)
imshow(im)
subplot(1,2,2)
imshow(rescaleIm(fidelity)*64,jet)

% Setup
nEvecs = min(nEvecs,nystromSize);
nystromSize = nystromSize + 1;

[E1, dex] = getMultimodalWeights(cat(3,im,lidar),[3 3],1,nystromSize,'L2');
[V1, D1] = getEvecs(E1, dex, nEvecs);

[E2, dex] = getTwoModalWeights(im,lidar,2,nystromSize,'L2',dex);
[V2, D2] = getEvecs(E2, dex, nEvecs);

[Einf, dex] = getTwoModalWeights(im,lidar,inf,nystromSize,'L2',dex);
[Vinf, Dinf] = getEvecs(Einf, dex, nEvecs);

[E0, dex] = getTwoModalWeights(im,lidar,0,nystromSize,'L2',dex);
[V0, D0] = getEvecs(E0, dex, nEvecs);

% V1 = matchEvecs(Vinf,V1);
% V2 = matchEvecs(Vinf,V2);
% V0 = matchEvecs(Vinf,V0);
% 
% %% Show the different V's
% for i = 1:5
%     subplot(2,2,1)
%     imshow(rescaleIm(reshape(V1(:,i),[size(im,1) size(im,2)]))*64,jet)
%     subplot(2,2,2)
%     imshow(rescaleIm(reshape(V2(:,i),[size(im,1) size(im,2)]))*64,jet)
%     subplot(2,2,3)
%     imshow(rescaleIm(reshape(Vinf(:,i),[size(im,1) size(im,2)]))*64,jet)
%     subplot(2,2,4)
%     imshow(rescaleIm(reshape(V0(:,i),[size(im,1) size(im,2)]))*64,jet)
%     pause
% end

% %% do spect clust on the V's
% s = rng('shuffle');
% K1 = eff_kmeans(V1,numClasses,200);
% rng(s)
% K2 = eff_kmeans(V2,numClasses,200);
% rng(s)
% Kinf = eff_kmeans(Vinf,numClasses,200);
% rng(s)
% K0 = eff_kmeans(V0,numClasses,200);
% 
% K1 = matchClasses(K1,Kinf);
% K2 = matchClasses(K2,Kinf);
% K0 = matchClasses(K0,Kinf);

%% Do MBO on the V's
K1 = MBO(V1, D1, numClasses, fidelity);
K2 = MBO(V2, D2, numClasses, fidelity);
Kinf = MBO(Vinf, Dinf, numClasses, fidelity);
K0 = MBO(V0, D0, numClasses, fidelity);

% imwrite(rescaleIm(reshape(K1,[size(im,1) size(im,2)]))*64,jet, ...
%     '~/Schoolwork/Chanussot/MeetingNotes/2018-01-24/1_norm.png')
% imwrite(rescaleIm(reshape(K2,[size(im,1) size(im,2)]))*64,jet, ...
%     '~/Schoolwork/Chanussot/MeetingNotes/2018-01-24/2_norm.png')
% imwrite(rescaleIm(reshape(Kinf,[size(im,1) size(im,2)]))*64,jet, ...
%     '~/Schoolwork/Chanussot/MeetingNotes/2018-01-24/max_norm.png')
% imwrite(rescaleIm(reshape(K0,[size(im,1) size(im,2)]))*64,jet, ...
%     '~/Schoolwork/Chanussot/MeetingNotes/2018-01-24/min_norm.png')


subplot(2,2,1)
imshow(rescaleIm(reshape(K1,[size(im,1) size(im,2)]))*64,jet)
subplot(2,2,2)
imshow(rescaleIm(reshape(K2,[size(im,1) size(im,2)]))*64,jet)
subplot(2,2,3)
imshow(rescaleIm(reshape(Kinf,[size(im,1) size(im,2)]))*64,jet)
subplot(2,2,4)
imshow(rescaleIm(reshape(K0,[size(im,1) size(im,2)]))*64,jet)

