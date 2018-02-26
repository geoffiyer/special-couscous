function [assign1to2, assign2to1] = MaxAssign(U, numMatches)
%MAXASSIGN Summary of this function goes here
%   Detailed explanation goes here

[~,tempU] = sort(U,2,'descend');
assign1to2 = [(1:size(tempU,1))', tempU(:,1:numMatches)];

[~,tempU] = sort(U,1,'descend');
assign2to1 = [(1:size(tempU,2))', tempU(1:numMatches,:)'];

end

