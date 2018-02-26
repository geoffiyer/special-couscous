%% Experimental graph match plan where we don't create an exact 1-1 match
% instead do something easier like a few top matches per node.
% Can stick this into the change detect code and I think it'll be okay??

function [assign1to2, assign2to1] = GraphMatch2(U1,U2,clusterSize)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%% Graph Laplacian and eigenvectors
% Eigenstuff using Nystrom

if(nargin < 3)
    clusterSize = 1;
end

%% Coarsify the graph

% Still don't have a good choice of how to coarsify.
[~,clusters1] = CoarsifyGraph(U1,clusterSize,'rand');
[~,clusters2] = CoarsifyGraph(U2,clusterSize,'rand');

%% First match the coarse points with each other

smallU = U1(clusters1(:,1),:)*(U2(clusters2(:,1),:)');
[smallAssign1to2,smallAssign2to1] = MaxAssign(smallU,4);

%% Now, for each coarse match (1 to 5), compare the corresponding patches to get the final match.

assign1to2 = CoarseAssign(U1,U2,clusters1,clusters2,smallAssign1to2,'max');
assign2to1 = CoarseAssign(U2,U1,clusters2,clusters1,smallAssign2to1,'max');

assign1to2 = sortrows(assign1to2);
assign2to1 = sortrows(assign2to1);
assign2to1 = [assign2to1(:,2), assign2to1(:,1)];

end