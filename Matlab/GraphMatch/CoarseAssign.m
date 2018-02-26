function [ assign ] = CoarseAssign( U1, U2, clusters1, clusters2, smallAssign, matchMethod)
%COARSEASSIGN Assign at coarse level first, then at fine level.
%   Input: U(i,j) = edgeweight between X1(i,:) and X2(j,:)
%          clusters1 = output from CoarsifyGraph(X1). Specifically, each
%                       row gives the indices of the corresponding cluster
%                       in X1 (so point i in Y1 corresponds to the many
%                       points clusters1(i,:) in X1)
%          clusters2 = similar for X2
%          smallAssign = assignments on the coarse level. Be careful of the indexing here.
%                        it is in terms of the small graph (so use
%                        clusters1, clusters2 to lift to big graph)
%          matchMethod = currently either 'hungarian' or 'max'
%
%   NOTE: Right now hungarian matching won't work unless coarse graph
%       cluster size exactly divides the full graph size.
%       (If not there will be some zeros in the clusters representing
%        unmatched pairs and I didn't write code for this yet).

if(~exist('matchMethod','var'))
    matchMethod = 'hungarian';
end


if(strcmp(matchMethod,'max'))    
    
    %% set up some variables
    assign = zeros(size(U1,1),2); % preallocate
    currentidx = 1;               % place tracker

    % For each match on the coarse data, extend to match on full data
    for i = 1:size(smallAssign,1)
        
        numMatches = size(smallAssign,2)-1;
        
        idx1 = clusters1(smallAssign(i,1),:)'; % indices of points in cluster i of X1
        idx2 = clusters2(smallAssign(i,2:(numMatches+1)),:); % indices of points in corresponding numMatches
        % worth of clusters from X2
        idx2 = idx2(:);
        
        idx1(idx1 == 0) = [];
        idx2(idx2 == 0) = [];
        
        Usmall = real(U1(idx1,:)*U2(idx2,:)'); % the corresponding section of U
        
        % Do the assign. Result is an indexing in terms of just this cluster
        % (not all of X1, X2). I call it local_assign.
        local_assign = MaxAssign(Usmall,1);
        
        % translate the indexing back to the full graph
        assign(currentidx:(currentidx+size(local_assign,1)-1),:) = ...
            [idx1(local_assign(:,1)), idx2(local_assign(:,2))];
        
        currentidx = currentidx + size(local_assign,1);
        
    end
else
    
    assign = zeros(min(size(U1,1),size(U2,1)),2);
    currentidx = 1;               % place tracker

    % For each match on the coarse data, extend to match on full data
    for i = 1:size(smallAssign,1)
        
        idx1 = clusters1(smallAssign(i,1),:)'; % indices of points in cluster i of X1
        idx2 = clusters2(smallAssign(i,2),:)'; % indices of points in cluster i of X2
        idx1(idx1 == 0) = [];
        idx2(idx2 == 0) = [];
        
        Usmall = real(U1(idx1,:)*U2(idx2,:)'); % the corresponding section of U
        
        % Do the assign. Result is an indexing in terms of just this cluster
        % (not all of X1, X2). I call it local_assign.
        local_assign = Hungarian(Usmall);
        
        % translate the indexing back to the full graph
        assign(currentidx:(currentidx+size(local_assign,1)-1),:) = ...
            [idx1(local_assign(:,1)), idx2(local_assign(:,2))];
        
        currentidx = currentidx + size(local_assign,1);
    end
end

end

