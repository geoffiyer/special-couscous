function [ E, dex ] = getMultimodalWeights( data, dims, p, m, normType, dex)
%GETWEIGHTS given the particular image/lidar inputs, give me the E
% weight matrix (nystrom style) that I'm going to use to get evecs.

addpath ../ImageTools

if(size(dims,1) > 1)
    dims = dims';
end
assert(size(dims,1) == 1,'what the heck why is dims a matrix??');
assert(sum(dims) == size(data,3),'Need dims to match data size');

data = reshape(data, [size(data,1)*size(data,2) size(data,3)]);
numSamples = size(data,1);
numModalities = size(dims,2);

if(nargin < 3)
    p = inf;
end
if(nargin < 4)
    m = min(101,size(data,1));
end
if(nargin < 5)
    normType = 'L2';
end
if(nargin <  6)
    % % Choose landmark points via a kmeans with few iterations
    % [~, center, m] = eff_kmeans(cat(2,data1,data2), m, 50); %#ite is restricted to 50
    % dex = dsearchn(cat(2,data1,data2),center);  % Find an actual pixel representative for each class
    
    % choose landmark points straight random
    dex = randperm(size(data,1));
end

centers = zeros(m, max(dims), numModalities);
dimidx = 1;
for i = 1:numModalities
    centers(:,1:dims(i),i) = data(dex(1:m),dimidx:(dimidx+dims(i)-1));
    dimidx = dimidx + dims(i);
end

%% W is already defined from E. Don't calculate twice (even though it doesn't matter)
% W = exp(-pNorm(sqrt(sqdist(center1', center1')),sqrt(sqdist(center2',center2')),p)/kernel.para);
Emat = zeros(numSamples,m,numModalities);
dimidx = 1;
for i=1:numModalities
    Emat(:,:,i) = sqrt(sqdist(data(:,dimidx:(dimidx+dims(i)-1))', centers(:,1:dims(i),i)'));
    dimidx = dimidx + dims(i);
end
scalings = zeros(1,numModalities);
for i=1:numModalities
    scalings(i) = sqrt(var(reshape(Emat(:,:,i),[numSamples*m 1])));
    Emat(:,:,i) = Emat(:,:,i)/scalings(i);
end
E = -pNorm(Emat,p);
E = exp(E/ (sqrt(var(E(:)))) );

end

