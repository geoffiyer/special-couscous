addpath ../Nystrom
addpath ../ImageTools

%% Initial input variables stuff whoopie
nystromSize = 100;
nEvecs = 50;
p = 2;

% %% The data
% numClasses = 3;
% im = [0; 2; 2.1; 2.2; 2.3; 2.2];
% lidar = zeros(size(im));
% fidelity = [1; 2; 2; 3; 3; 0];
% fidelity = fidelity/255;
% fidelityConst = 1.0e0;
% dtConst = 1.0e0;

numClasses = 6;
A = load('../SpectralClustering/Umbrella/Data.mat');
im = A.im;
lidar = A.lidar;
lidar = rescaleIm(lidar);
im = rescaleIm(im);
fidelity = double(imread('~/Schoolwork/Chanussot/MBO_Code/data/fidelity_umbrella_6class.tiff'));
fidelity = fidelity/255; % rescale to handle how matlab prints image to text
fidelityConst = 1.0e1;
dtConst = 1.0;

% numClasses = 4;
% im = zeros(150,150,3);
% lidar = zeros(150,150,3);
% 
% for i = 1:size(im,1)
%     for j = 1:size(im,2)
%         if i <= size(im,1)/2
%             im(i,j,:) = [1 1 0];
%         else
%             im(i,j,:) = [0 1 1];
%         end
%         if j <= size(im,2)/2
%             lidar(i,j,:) = [1 1 0];
%         else
%             lidar(i,j,:) = [0 1 1];
%         end
%     end
% end
% im = im + randn(size(im,1),size(im,2))*0.002;
% lidar = lidar + randn(size(im,1),size(im,2))*0.002;
% 
% fidelity = zeros([size(im,1) size(im,2)]);
% fidelity(40:45,40:45) = 1/255;
% fidelity(100:105,100:105) = 2/255;
% fidelity(40:45,100:105) = 3/255;
% fidelity(100:105,40:45) = 4/255;
% 
% subplot(1,2,1)
% imshow(im)
% subplot(1,2,2)
% imshow(rescaleIm(fidelity)*64,jet)

%% Setup
nEvecs = min(nEvecs,nystromSize);
nystromSize = nystromSize + 1;

[E, dex] = getTwoModalWeights(im,lidar,2,nystromSize,'L2');
[V, D] = getEvecs(E, dex, nEvecs);
V = orth(V);

% [Einf, dex] = getTwoModalWeights(im,lidar,inf,nystromSize,'L2',dex);
% [Vinf, Dinf] = getEvecs(Einf, dex, nEvecs);

%% Do MBO on the V's
[K, info] = MBO(V, D, numClasses, fidelity, fidelityConst, dtConst);

if(false)
% imwrite(rescaleIm(reshape(K1,[size(im,1) size(im,2)]))*64,jet, ...
%     '~/Schoolwork/Chanussot/MeetingNotes/2018-01-24/1_norm.png')
% imwrite(rescaleIm(reshape(K2,[size(im,1) size(im,2)]))*64,jet, ...
%     '~/Schoolwork/Chanussot/MeetingNotes/2018-01-24/2_norm.png')
% imwrite(rescaleIm(reshape(Kinf,[size(im,1) size(im,2)]))*64,jet, ...
%     '~/Schoolwork/Chanussot/MeetingNotes/2018-01-24/max_norm.png')
% imwrite(rescaleIm(reshape(K0,[size(im,1) size(im,2)]))*64,jet, ...
%     '~/Schoolwork/Chanussot/MeetingNotes/2018-01-24/min_norm.png')
end

subplot(1,2,1)
imshow(reshape((K+1)/(numClasses),[size(im,1) size(im,2)])*64,jet)
subplot(1,2,2)
jetshow(fidelity)
