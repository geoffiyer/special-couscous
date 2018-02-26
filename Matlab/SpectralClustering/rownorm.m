function [ B ] = rownorm( A )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

B = sqrt(sum(A.^2,2));

end

