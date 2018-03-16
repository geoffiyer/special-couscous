function [ V, va ] = getEvecs( E, dex, nEvecs )
%GETEVECS Use the nystrom-style weight matrix to get the evecs
%   So the weights have already been created.
%   This function just does the linear algebra part,
%   Then returns some evecs.
%   INPUT: E = the (# points) x (# landmark) weight matrix
%          dex = the indeces of landmark points
%          nEvecs = number of evecs to grab (optional parameter)
%                   can't be bigger than size(E,2)
%   OUTPUT: V = the (# points) x (# evecs) eigenvector matrix

%% The actual linear algebra
m = size(E,2);
n = size(E,1);
W = E(dex(1:m),:);  % m x m matrix

G = E * W^(-1/2); % m x n * m x m -> m x n

d = G * (G' * ones(n,1));
% Geoff says: I changed the below line. Original is d(find(d<0)) = 1e-5;
d(d < 0) = 1e-15;
d = d.^(-0.5);
%%%%%%%%%%%%%%%%%%%
% Old code. changes each row of G by left-multiplying by a diagonal mat
%  I think it's smarter to just scale each row, as you can see below.
% sd = sparse(n, n);
% for i = 1:n
%     sd(i,i) = d(i);
% end;
% G = sd * G;
%%%%%%%%%%%%%%%%%%%
for i = 1:n
    G(i,:) = d(i) * G(i,:);
end
[U, L] = eig(G'*G);
va = diag(L);
[~, dex] = sort(va,'descend');
U = U(:,dex);
V = G * U;
for i = 1:n  % this loop replaces  V = sd*V
    V(i,:) = d(i) * V(i,:);
end

V = real(V);
% V = U;

%% Figure out how many evecs to use
% Do a 2-class threshold on the evals D, use this to pick the cutoff

if(nargin < 3)
    evecsToUse = 2-otsu(va,2);
    if(evecsToUse(1) == 0)
        evecsToUse = 1-evecsToUse;
    end
    V = V(:,(evecsToUse == 1));
    V = V(:,2:size(V,2));
    va = va(1:size(V,2));
else
    nEvecs = min(nEvecs, size(V,2)-1);
    V = V(:,2:(nEvecs+1));
    va = va(1:size(V,2));
end

end
