function [ E, dex ] = getMultimodalWeights( data, dims, p, m, normType, dex)
%GETWEIGHTS given the particular image/lidar inputs, give me the E
% weight matrix (nystrom style) that I'm going to use to get evecs.

% the normType variable is gonna be confusing
% within each modality, there are two choices for norm. 2-norm and angle.
% normType = 1 x numModalities vector
%  1 at position i means use angle norm for modality i.
%  0 means use 2 norm
% Another note: A big reason for why is speed. These square dist matrices
% are pretty big! We want cool matlab-y tricks to compute.

addpath ../ImageTools

if(size(data,3) > 1)
    warning('I changed the code you should reshape data to 2D-array before getting weights');
end
if(size(dims,1) > 1)
    dims = dims';
end
assert(size(dims,1) == 1,'what the heck why is dims a matrix??');
assert(sum(dims) == size(data,2),'Need dims to match data size');

numSamples = size(data,1);
numModalities = size(dims,2);

if(nargin < 3)
    p = inf;
end
if(nargin < 4)
    m = min(101,size(data,1));
end
if(nargin < 5)
    % If no input, then all intra-modality norms are just 2-norm
    normType = zeros(size(dims));
end
if(nargin <  6)
    % % Choose landmark points via a kmeans with few iterations
    % [~, center, m] = eff_kmeans(cat(2,data1,data2), m, 50); %#ite is restricted to 50
    % dex = dsearchn(cat(2,data1,data2),center);  % Find an actual pixel representative for each class
    
    % choose landmark points straight random
    dex = randperm(size(data,1));
end

assert(isequal(size(dims),size(normType)),'need one normType for each modality');

centers = zeros(m, max(dims), numModalities);
dimidx = 1;
for i = 1:numModalities
    centers(:,1:dims(i),i) = data(dex(1:m),dimidx:(dimidx+dims(i)-1));
    dimidx = dimidx + dims(i);
end

% We're getting a real problem with space complexity. Storing 8 different
% copies of E (one for each modality) is too hard. Instead, store 2. One is
% the final, the other is the dude we're working on right now. Each time we
% make a new dude we update the final.
E = zeros(numSamples,m);
scalings = zeros(1,numModalities);
dimidx = 1;
for i=1:numModalities
    if(normType(i) == 0)
        Enew = sqrt(sqdist(data(:,dimidx:(dimidx+dims(i)-1))', centers(:,1:dims(i),i)'));
    else
        Enew = angledist(data(:,dimidx:(dimidx+dims(i)-1)), centers(:,1:dims(i),i));
    end
    dimidx = dimidx + dims(i);
    scalings(i) = sqrt(var(reshape(Enew,[numSamples*m 1])));
    Enew = Enew/scalings(i);
    % E_update = E_old + E+new^p
    if p == inf
        E = max(E,Enew);
    elseif p == 0
        E = min(E,Enew);
    else
        E = E + Enew.^p;
    end
end
    
if p ~= inf && p ~= 0
    E = -(E.^(1/p));
end
E = exp(E/ (sqrt(var(E(:)))) );

end