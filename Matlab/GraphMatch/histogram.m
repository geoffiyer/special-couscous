function [ H1,H2 ] = histogram(U1, U2, B)
% Okay, U1 and U2 are matrices of eigenvectors (each column is 1 evec)
% B is the number of boxes for this histogram
% I don't know why I need the same number of evecs in each U that should
% probably go.
% 

assert(size(U1,2) == size(U2,2), 'Need same number of evecs in each U');
k = size(U1,2);

if (nargin < 3)
    B = 100;
end

globalmin = min(min([U1;U2]));
globalmax = max(max([U1;U2]));
delta = (globalmax-globalmin)/B;

H1 = zeros(B,k);
H2 = zeros(B,k);

for i=1:B
    low = globalmin + (i-1)*delta;
    high = low + delta;
    H1(i,:) = sum( U1 >= low & U1 < high,1);
    H2(i,:) = sum( U2 >= low & U2 < high,1); 
end

end