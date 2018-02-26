function [ W ] = weightMatrix( data, p, sigma )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if(nargin < 2)
    p = 2;
end

n = size(data,1);
W = zeros(n);

if(p == Inf)
    for i = 1:n-1
        for j = i:n
            W(i,j) = max(abs(data(i,:) - data(j,:)));
        end
    end
elseif(p == 0)
    for i = 1:n-1
        for j = i:n
            W(i,j) = min(abs(data(i,:) - data(j,:)));
        end
    end
else
    for i = 1:n-1
        for j = i:n
            temp = abs(data(i,:) - data(j,:));
            W(i,j) = sum(temp.^p)^(1/p);
        end
    end
end

if(nargin < 3)
    sigma = std(W(W ~= 0));
    sigma = 1;
end

W = exp((W + W')/-sigma);

end