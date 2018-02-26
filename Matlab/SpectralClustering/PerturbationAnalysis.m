
addpath ../ImageTools

% A = load('DFC/DFC2013.mat');
% baseIm = A.im(:,:,[40,70,100]);
% nclasses = 6;

% A = load('Jadeplant/Data.mat');
% unionClasses1 = 3;
% unionClasses2 = 2;
% nclasses = 7;

A = load('Umbrella/Data.mat');
nclasses = 6;
epsilon = 0.025;

% A = load('DFC/DFC2015.mat'); % For DFC, nclasses = 5 looks pretty good.
%                              % Should probably think more about how to
%                              % choose number.
% nclasses = 6;

p = inf;
numEvecs = 30;
numTrials = 30;

%% NOTE TO SELF FROM THIS TEST:
% Okay, so the final K-means for spectral clustering has a HUGE variation
% This really sucks to deal with
% So I added a bit to make the random seed the same for both tests
% for just this last kmeans part
% Remaining elements:
%  kmeans to choose landmark nodes for Nystrom: tested and OKAY
%  adding an epsilon rng to weight matrix: not yet tested we'll see
%      Specifically we should start looking at I (the evecs) after an epsilon

KpercentList = zeros(numTrials,1);
VList = zeros(numTrials,numEvecs);
[E,dex] = getWeights(A.im, A.lidar, p);
V = getEvecs(E,dex);
V = V./colnorm(V);

tic
for i = 1:numTrials
%     rng('shuffle')
%     s = rng;
%     [error,I,K,~,~,E] = SpectralClassification(A.im, A.lidar, nclasses, displayImage, smoothing, p, 0, s);
%     [errorPert,Ipert,Kpert,~,~,Epert] = SpectralClassification(A.im, A.lidar, nclasses, displayImage, smoothing, p, epsilon, s);
% 
%     Kpert = matchClasses(Kpert,K);
%     KchangePercent = sum(sum(K ~= Kpert))/(size(K,1) * size(K,2));
%     KpercentList(i) = KchangePercent;
    Epert = E + randn(size(E))*std(E(:))*epsilon;
    Vpert = getEvecs(Epert,dex,size(V,2));
    Vpert = Vpert./colnorm(Vpert);
    Vpert = matchEvecs(V,Vpert);
    nevecsToUse = min([numEvecs,size(V,2),size(Vpert,2)]);
    VList(i,1:nevecsToUse) = reshape(min( colnorm(V(:,1:nevecsToUse) - Vpert(:,1:nevecsToUse)), ...
                                          colnorm(V(:,1:nevecsToUse) + Vpert(:,1:nevecsToUse))) ...
                                       ./ colnorm(V(:,1:nevecsToUse)) , [1 nevecsToUse]);
    if(VList(i,1) > 0.5)
        Vsave = Vpert;
        idxSave = i;
    end
end
runtime = toc/numTrials;
imagesc(VList)

% I = reshape(V,[size(A.im,1) size(A.im,2) size(V,2)]);
% Ipert = reshape(Vsave,[size(A.im,1) size(A.im,2) size(Vsave,2)]);
% for i = 1:nevecsToUse
%     subplot(1,2,1);
%     imshow(rescaleIm(I(:,:,i))*64,jet);
%     subplot(1,2,2);
%     imshow(rescaleIm(Ipert(:,:,i))*64,jet);
%     pause
% end
% clear i

% subplot(1,3,1)
% imshow(rescaleIm(K)*64,jet)
% title('True Weights')
% subplot(1,3,2)
% imshow(rescaleIm(Kpert)*64,jet)
% title('Perturbed Weights')
% subplot(1,3,3)
% imshow((K ~= Kpert))
% title('changed pixels')