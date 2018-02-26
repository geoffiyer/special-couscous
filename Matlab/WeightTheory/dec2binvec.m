function [ vec ] = dec2binvec( n, length )
%DEC2BIN Summary of this function goes here
%   Detailed explanation goes here

vec = zeros(floor(log2(n))+1,1);

for i = 1:size(vec)
    vec(i) = mod(n,2);
    n = floor(n/2);
end

if(nargin == 2)
    vec = [vec; zeros(length - size(vec,1),1)];
end

vec = logical(vec);

end

