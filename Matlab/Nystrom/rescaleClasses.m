function [ Knew ] = rescaleClasses( K, killLevel )
%RESCALEFIDELITY Summary of this function goes here
%   Detailed explanation goes here

if(nargin < 2)
    killLevel = inf;
end

nums = [unique(K), zeros(size(unique(K),1),1)];

for i = 1:size(nums,1)
    nums(i,2) = size(K(K == nums(i,1)),1);
end
nums = sortrows(nums,-2);
nums(nums(:,2) < size(K(:),1)/killLevel, :) = [];

Knew = zeros(size(K));

for i = 1:size(nums)
    Knew(K == nums(i,1)) = i - 1;
end

end

