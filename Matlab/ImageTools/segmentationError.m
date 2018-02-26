function [ e ] = segmentationError( data, K, mode)
% For each class (resp. region), compare actual data against average value
% over the class (resp. region). Add up the norm of the difference to find
% the segmentation error.
%   Data = data
%   K = classification
%   mode = either 'classes' or 'regions'

% If they didn't input a mode, do it via classes
if(~exist('mode','var'))
    mode = 'classes';
end

if(strcmp(mode,'regions'))
    K = connComp(K);
end

data = reshape(data, [size(data,1)*size(data,2) size(data,3)]);
K = reshape(K, [size(K,1)*size(K,2) 1]);

nclasses = max(max(K));

e = 0;

for i = 1:nclasses
    idx = (K == i);
    idxsize = sum(idx);
    diffmat = data(idx,:) - repmat(mean(data(idx,:)), [idxsize 1]);
    e = e + sum(sqrt(sum(diffmat.^2,2)));
end

e = e / ( size(K,1)*size(K,2));

end