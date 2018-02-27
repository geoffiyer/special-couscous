function [ E, dex ] = getTwoModalWeights( im, lidar, p, m, normType, dex)
%GETWEIGHTS given the particular image/lidar inputs, give me the E
% weight matrix (nystrom style) that I'm going to use to get evecs.

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

if(nargin < 3)
    p = inf;
end
if(nargin < 4)
    m = 101;
end
if(nargin < 5)
    normType = 'L2';
end
if(nargin <  6)
    % % Choose landmark points via a kmeans with few iterations
    % [~, center, m] = eff_kmeans(cat(2,data1,data2), m, 50); %#ite is restricted to 50
    % dex = dsearchn(cat(2,data1,data2),center);  % Find an actual pixel representative for each class
    
    % choose landmark points straight random
    dex = randperm(size(X1,1));
end

m = min(m, size(X1,1));

center1 = X1(dex(1:m),:);
center2 = X2(dex(1:m),:);

%% W is already defined from E. Don't calculate twice (even though it doesn't matter)
% W = exp(-pNorm(sqrt(sqdist(center1', center1')),sqrt(sqdist(center2',center2')),p)/kernel.para);
if(strcmp(normType,'angle'))
    % data1 = nonlocal(data1,4,3,imSize);  % Some nonlocal thing, but it's
    % sketchily written...
    % This is likely trash.
    E1 = angledist(X1,center1);
else
    E1 = sqrt(sqdist(X1', center1'));
end
E2 = sqrt(sqdist(X2', center2'));
scaling1 = sqrt(var(E1(:)));
scaling2 = sqrt(var(E2(:)));
E1 = E1 / scaling1;
E2 = E2 / scaling2;
% Nonlocal means part: currently not working
% squaresize = 1;
% E1 = nonlocal(E1, squaresize, imSize);
% E2 = nonlocal(E2, squaresize, imSize);
E = -pNorm(cat(3,E1,E2),p);
E = exp(E/ (sqrt(var(E(:)))) );

end

