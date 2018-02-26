function [ W ] = WeightMatrix( X, type_string, numNeighbors )
% Calculate weight matrix using either L2 norm or knn
%   X = data matrix. #rows = #data points. #cols = dim of data
%   type_string = either 'L2' or 'knn'
%   numNeighbors = which k to use for knn. Only used if type is knn

if(nargin < 2)
    type_string = 'L2';
elseif (strcmp(type_string,'L2') && nargin < 3)
    numNeighbors = 20;
end

if (strcmp(type_string,'knn'))
    neighborIds = kNearestNeighbors(X,X,numNeighbors);
    W = 3*ones(size(X,1));
    for i = 1:size(neighborIds,1)
        W(i,neighborIds(i,:)) = 0;
    end
    W = (W + W')/2;
    W = exp(-W);
else
    W = sqrt(sqdist(X',X'));
    W = W/std(W(:));
    W = exp(-W);
end

end

