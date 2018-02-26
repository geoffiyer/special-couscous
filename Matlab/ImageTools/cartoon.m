function [ result ] = cartoon( data, K, mode )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if ~exist('mode','var')
    mode = 'regions';
end

if (strcmp(mode,'regions'))
    labels = connComp(K);
elseif (strcmp(mode,'classes'))
    labels = K;
else
    disp('Nonsense input to cartoon function');
    return;
end

originalSize = size(data);

data = reshape(data, [size(data,1)*size(data,2) size(data,3)]);
labels = reshape(labels, [size(labels,1)*size(labels,2) 1]);
result = zeros(size(data));

for i = 1:(max(max(labels)))
    idx = (labels == i);
    count = sum(idx);
    result(idx,:) = repmat(mean(data(idx,:)),[count 1]);
end

result = reshape(result,originalSize);

end