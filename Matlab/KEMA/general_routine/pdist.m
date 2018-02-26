function [ D ] = pdist( X )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

m = size(X,1);
D = zeros(1,m*(m-1)/2);
count = 1;

for i=1:m
    for j = (i+1):m
        D(count) = norm(X(j,:) - X(i,:));
        count = count+1;
    end
end

end

