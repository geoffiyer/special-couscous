function [ A ] = pNorm( data, p )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

if p == inf
    A = max(data,[],3);
elseif p == 0
    A = min(data,[],3);
else
    A = sum(data.^p,3).^(1/p);
end
end

