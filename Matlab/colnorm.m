function [ B ] = colnorm( A )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

B = sqrt(sum(A.^2,1));

end

