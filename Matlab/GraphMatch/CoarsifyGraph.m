function [Y, clusters] = CoarsifyGraph(X,numNeighbors,method)
%% Take a dataset X represented by a weighted graph (matrix W).
% Find a subgraph of X that's supposed to be a coarse representation
% Input: X, a (num_points) x (point_dim) matrix
%        numNeighbors, the amount by which to scale down the graph
% Output: Y, the coarse version of X
%         clusters, a matrix keeping track of how indices in Y
%             match in X
%            Row i of clusters will contain the many indices in X
%            that match with index i in Y

if( nargin < 2)
    numNeighbors = 10;
end

if( nargin < 3)
    method = 'degree';
end



if( strcmp(method,'degree'))
    
    % Calculate degrees
    W = sqrt(sqdist(X',X'));
    W = W/std(W(:));
    W = exp(-W);
    D = sum(W,1);
    
    
    %% Build Y iteratively
    % Grab the most influential node from X. Remove it and it's
    % (numNeighbors-1) nearest neighbors. Repeat.
    Y = zeros(floor(size(X,1)/numNeighbors), size(X,2));
    clusters = zeros(size(Y,1),numNeighbors);
    keep = 1:size(X,1);
    count = 1;
    
    while( size(keep,2) >= numNeighbors)
        [~,maxInd] = max(D);
        neighbors = kNearestNeighbors(X(keep,:),X(maxInd,:),numNeighbors);
        Y(count,:) = X(maxInd,:);
        clusters(count,:) = keep(neighbors);
        D(keep(neighbors)) = -inf;
        keep(neighbors) = [];
        count = count + 1;
    end
    
    if(size(keep,2) ~= 0)
        clusters(count,:) = [keep, zeros(1,numNeighbors - size(keep,2))];
    end
    
elseif(strcmp(method,'rand'))
    
    %% Build Y iteratively
    % Pick a random node from X. Remove it and it's
    % (numNeighbors-1) nearest neighbors. Repeat.
    Y = zeros(ceil(size(X,1)/numNeighbors), size(X,2));
    clusters = zeros(size(Y,1),numNeighbors);
    keep = 1:size(X,1);
    count = 1;
    
    while( size(keep,2) >= numNeighbors)
        ind = keep(randi(size(keep,2)));
        neighbors = kNearestNeighbors(X(keep,:),X(ind,:),numNeighbors);
        Y(count,:) = X(ind,:);
        clusters(count,:) = keep(neighbors);
        keep(neighbors) = [];
        count = count + 1;
    end
    
    if(size(keep,2) ~= 0)
        clusters(count,:) = [keep, zeros(1,numNeighbors - size(keep,2))];
        Y(count,:) = X(keep(1),:);
    end
    
else
    %% Use kmeans to pick out the right number of centers
    % These are the same as the 'influential nodes' above
    % For each center, make a cluster out of its (numNeighbors) nearest
    % neighbors.
    
    addpath('./Nystrom/');
    [idx,Y,m] = eff_kmeans(X, round(size(X,1)/numNeighbors), 100);
    clusters = cell(m,1);
    for i = 1:m
        clusters{i} = find(idx == i);
    end
    rmpath('./Nystrom/');
    
end

end