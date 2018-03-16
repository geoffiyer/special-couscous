function [ GTnew, idx ] = rescaleFidelity( GT, killLevel )
% If GTnew(i,j) == k, then GT(i,j) = idx(k)
%   If you think of idx as a function {1:n} -> {1:n}, then
%   idx(GTnew) = GT. Of course matlab syntax doesn't work that way

if(nargin < 2)
    killLevel = inf;
end

nums = [unique(GT), zeros(size(unique(GT),1),1)];
nums(nums(:,1) == 0,:) = [];

for i = 1:size(nums,1)
    nums(i,2) = size(GT(GT == nums(i,1)),1);
end
[nums, sortidx] = sortrows(nums,-2);  % nums_new = nums_old(idx,:)
sortidx(nums(:,2) < size(GT(:),1)/killLevel, :) = [];
nums(nums(:,2) < size(GT(:),1)/killLevel, :) = [];

GTnew = zeros(size(GT));
idx = zeros(size(nums,1),1);

for i = 1:size(nums)
    GTnew(GT == nums(i,1)) = i;
    idx(i) = nums(i,1);
end

end

