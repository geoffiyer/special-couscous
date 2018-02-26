%% This code implements the example we talked about in our meeting 11-15.

% The ground truth data consists of two clouds. The main cloud is centered
% around a line, and the secondary cloud is a cluster.
% In the dataset X1, the two clouds are merged. In X2 they are separated.

%% Create the dataset
n1 = 40;    % size of main cloud
n2 = 10;     % size of secondary cloud

xmin = 0;   % coordinates of line
xmax = 20;
ymean = 0;
ysdev = 1;  % distribution of y-values around line

% Data between X1 and X2 is assumed to be co-registered. But we add a small
% Gaussian noise to simulate differences in observations between datasets.
noisesdev = 1;

% distribution of secondary cloud
cloud2xmean = 13;
cloud2xsdev = 1;
cloud2ymean = 7;
cloud2ysdev = ysdev;

% initialize data sets to zero
X1 = zeros(2,n1+n2);
X2 = zeros(2,n1+n2);

% Create x-values for both clouds. The points on the line are random
% uniform. The points in the cluster are normally distributed.
X1(1,1:n1) = rand(1,n1)*(xmax-xmin) + xmin;
X1(1,(1+n1):(n1+n2)) = randn(1,n2)*cloud2xsdev+cloud2xmean;
X2(1,:) = X1(1,:) + randn(1,n1+n2)*noisesdev;

% Create y-values for both clouds. In X1 the y-values all come from the
% same distribution. In X2 the points from the secondary cloud have a
% different mean.
X1(2,:) = randn(1,n1+n2)*ysdev + ymean;
X2(2,:) = X1(2,:) + randn(1,n1+n2)*noisesdev;
X2(2,(1+n1):(n1+n2)) = X2(2,(1+n1):(n1+n2)) + (cloud2ymean - ymean);

% Label our clouds.
cloud1 = (1:n1);
cloud2 = (n1+1):(n1+n2);

%% Create the Graph-Laplacian. Find Eigenvalues and Eigenvectors
% Note: This code isn't optimized.
[nrx1, ncx1] = size(X1);
[nrx2, ncx2] = size(X2);
W = zeros(ncx1,ncx1);   % Weight matrix
for i = 1:ncx1
    for j = (i+1):ncx1
        W(i,j) = exp(-1*max( norm(X1(:,i)-X1(:,j)), norm(X2(:,i)-X2(:,j))));
        W(j,i) = W(i,j);
    end
end

L = diag(sum(W,2)) - W; % Graph Laplacian
[V,D] = eig(L);     % Eigenvalues/Eigenvectors of Graph Laplacian

%% Plot the data
figure
hold on;
scatter(X1(1,cloud1),X1(2,cloud1), '+', 'red');
scatter(X1(1,cloud2),X1(2,cloud2), 'o', 'blue');
title('X1 Dataset');
legend('Cloud 1','Cloud 2','Location','bestoutside');
hold off;

figure
hold on;
scatter(X2(1,cloud1),X2(2,cloud1), '+', 'red');
scatter(X2(1,cloud2),X2(2,cloud2), 'o', 'blue');
title('X2 Dataset');
legend('Cloud 1','Cloud 2','Location','bestoutside');
hold off;

figure
hold on;
scatter(V(cloud1,2),V(cloud1,3), '+', 'red');
scatter(V(cloud2,2),V(cloud2,3), 'o', 'blue');
title('Latent Space (using 2 eigenvectors)');
legend('Cloud 1','Cloud 2','Location','bestoutside');
hold off;