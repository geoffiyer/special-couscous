function [error,I,K,scaling1,scaling2,E] = SpectralClassification(im, lidar, nclasses, displayImage, smoothing, p, epsilon, seed)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

addpath ../ImageTools

nrow = size(im,1);
ncol = size(im,2);
if (nrow ~= size(lidar,1) || ncol ~= size(lidar,2))
    disp('image and lidar size do not match');
    return
end

X1dim = size(im,3);
X2dim = size(lidar,3);
X1 = reshape(im,[nrow*ncol X1dim]);
X2 = reshape(lidar, [nrow*ncol X2dim]);

if ~exist('seed','var')
    rng('shuffle')
end
if ~exist('epsilon','var')
    epsilon = 0;
end
if ~exist('p','var')
    p = inf;
end
if ~exist('smoothing','var')
    smoothing = 1;
end
if ~exist('displayImage','var')
    displayImage = false;
end

%m = floor( sqrt(sqrt(nrow*ncol)) *3 ) ;  % number of landmark points. m << Number of Data Points
m = min(101,size(X1,1));

% Right now 'para' doesn't do anything. It's supposed to have to do with
% the scaling of the distance matrix. Maybe do that later.
[V,D,scaling1,scaling2,E] = INys_SpectrEmbed(X1,X2,m, p,'L2',epsilon);         % V = first m eigenvectors of the graph laplacian

%% Display some results.

I = reshape(V(:,2:size(D,1)),[nrow ncol (size(D,1)-1)]);
%I = (I - min(min(I)))./ (max(max(I)) - min(min(I)));

%% Do some kmeans classification and display those results

% Old method: Use all eigenvectors
% nevecsToUse = m-1;

% assignin('base','D',diag(D)); %temporary output for debugging

% New method: Use until evals get to big (haven't optimized any constants
% yet)
b = max(max(D(1:3,1:3)));
c = find(diag(D)>b*4,1);
nevecsToUse = c-1;
assignin('base','nevecsUsed',nevecsToUse);

if exist('seed','var')
    rng(seed)
end
maxIter = 150;
[idx, ~, ~ ] = eff_kmeans(V(:,2:(1+nevecsToUse)),nclasses,maxIter);
K = reshape(idx, [nrow, ncol]);
K = smoothPixels(K,smoothing);

% error = segmentationError(cat(3,im/scaling1,lidar/scaling2),K);  %
% segmentation by classes

error = segmentationError(cat(3,im/scaling1,lidar/scaling2),K);

if displayImage
    mainfig = figure;
    mainfig.Position = [300,100,1200,850];
    if(size(im,3) == 3)
        subplot(2,2,1)
        imshow(im);
        title('Original Optical Image');
    end
    subplot(2,2,2);
    if(size(lidar,3) > 1)
        imshow(lidar/max(max(max(lidar))));
    else
        imshow(lidar/max(max(lidar)));
    end
    title('Original Lidar Image');
    subplot(2,2,3);
    imagesc(K);
    title('Kmeans Classification');
    subplot(2,2,4);
    imshow(imWithBorders);
    title('Classification on Optical Image');
end

%for i = 1:nevecsToDisplay
%    figure
%    imshow(I(:,:,i));
%    title(strcat('Eigenvector number ',num2str(i)))
%end

end