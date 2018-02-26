function [ U1, U2 ] = getEigenvectors( X1, X2, normType )
%GETEIGENVECTORS Summary of this function goes here
%   Detailed explanation goes here

addpath('../Nystrom/');
n = max(size(X1,1),size(X2,1));
m = max(round(sqrt(sqrt(n))*4));
m = max(100,m);
m = min([size(X1,1),size(X2,1),m]);
[U1,D1] = INys_SpectrEmbed(X1,m,normType);
[U2,D2] = INys_SpectrEmbed(X2,m,normType);
U1 = real(U1);
U2 = real(U2);
rmpath('../Nystrom/');

% How many eigenvectors to use?
% I decided on using all between smallest and 2*smallest
% this hasn't been optimized though. It's pretty much done at random
b1 = min(max(D1(2:size(D1,1),2:size(D1,2))));
c1 = find(diag(D1)>b1*2,1)-2;
b2 = min(max(D2(2:size(D2,1),2:size(D2,2))));
c2 = find(diag(D2)>b2*2,1)-2;
nevecsToUse = floor(max([c1,c2]));
nevecsToUse = min(n,nevecsToUse);
nevecsToUse = max(2,nevecsToUse);

U1 = U1(:,2:(2+nevecsToUse-1));
U2 = U2(:,2:(2+nevecsToUse-1));

% Try normalizing the rows.
U1 = U1 ./ sqrt(sum(U1.^2,2));
U2 = U2 ./ sqrt(sum(U2.^2,2));


end

