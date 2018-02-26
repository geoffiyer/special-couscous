function [ A ] = pNorm( X1, X2, p )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

if p == inf
    A = max(X1,X2);
elseif p == -inf
    A = min(X1,X2);
else
    A = (X1.^p + X2.^p).^(1/p);
end
end

