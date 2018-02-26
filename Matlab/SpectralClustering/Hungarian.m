function [ assign ] = Hungarian( U )
%HUNGARIAN Uses hungarian algorithm to determine an matching that minimizes
%total cost.
%   assign will be an n x 2 matrix.
%      n = min dimension of U
%   Each row gives a match between set1 and set2

[rowAssign,~] = munkres(1-U); 
assign = [(1:size(rowAssign,2))' , rowAssign'];
zero_rows = (min(assign,[],2) == 0);
assign(zero_rows,:) = [];

end

