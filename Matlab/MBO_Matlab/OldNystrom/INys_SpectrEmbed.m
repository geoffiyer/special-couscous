% Improved Nystrom method for Spectral Embedding
% This function implements the improved Nystrom method in <Improved Nystrom
% Low Rank Approximation and Error Analysis>, which can be used for Laplacian 
% Eigenmap, Spectral Clustering, or Normalized cut.

% Input:
% data1: n-by-dim1 data matrix. Co-registered with data2;
% data2: n-by-dim2 data matrix. Co-registered with data1;
% kernel: a struct with two elements;
%         kernel.type: 'pol' or 'rbf';
%         kernel.para: d in the polynomial kernel <x,y>^d;
%                      b in the rbf kernel exp(-||x||^2/b);
% m: number of landmark points, ucually chosen much smaller than data size n.

% Output:
% V: n-by-m matrix, containing the top m eigenvectors of the normalized kernel matrix 
% D^(-0.5) * K * D^(-0.5), and these eigenvectors are further normalized by D^{-1/2} as
% in the normalized cut method; here K is the kernel matrix and D is the degree matrix. 
% The eigenvectors are sorted in descending order by corresponding eigenvalues.
% *note* you might want to re-scale each column of V for further processing

function V = INys_SpectrEmbed(data1, data2, kernel, m)

[n, dim1] = size(data1);
dim2 = size(data2,2);

%% Geoff says: This code likes to choose landmark points with kmeans.
% I don't know a smart way to do this with 2 data sets so I'm just going
% to choose landmark points randomly (the Nystrom I'm used to does this).

% Their kmeans code:
% [idx, center, m] = eff_kmeans(data, m, 5); %#ite is restricted to 5

%% random sampling of landmark points
dex = randperm(n);
center1 = data1(dex(1:m),:);
center2 = data2(dex(1:m),:);

if(kernel.type == 'pol')
    W = center * center';
    E = data * center';
    W = W.^kernel.para;
    E = E.^kernel.para;
end

% Geoff says: My code
if(kernel.type == 'rbf')
    W = exp(-max(sqrt(sqdist(center1', center1')),sqrt(sqdist(center2',center2')))/kernel.para);
    E = exp(-max(sqrt(sqdist(data1', center1')),sqrt(sqdist(data2',center2')))/kernel.para);
end;

% Geoff says: Their original code
% if(kernel.type == 'rbf')
%     W = exp(-sqdist(center1', center1')/kernel.para);
%     E = exp(-sqdist(data1', center1')/kernel.para);
% end;

G = E * W^(-1/2);

d = G * (G' * ones(n,1));
% Geoff says: I changed the below line. Original is d(find(d<0)) = 1e-5;
d(d < 0) = 1e-15;
d = d.^(-0.5);
% Geoff says: The part with sd can be improved for sure. Why build a sparse
% diagonal matrix to left-multiply when I can just multiply each column?
sd = sparse(n, n);
for i = 1:n
    sd(i,i) = d(i);
end;
G = sd * G;
[U, L] = eig(G'*G);
va = diag(L);
[va, dex] = sort(va,'descend');
U = U(:,dex);
V = G * U;
V = sd * V;
