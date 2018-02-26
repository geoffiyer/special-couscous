
%% This will print the thing I've been calling set1 inside the generatedata thingy

addpath ../Matlab/

% Print ground truth
figure
hold on;
scatter3(XG(cloud1,1),XG(cloud1,2),XG(cloud1,3),'o','blue');
scatter3(XG(cloud2,1),XG(cloud2,2),XG(cloud2,3),'+','red');
scatter3(XG(cloud3,1),XG(cloud3,2),XG(cloud3,3),'x','black');
legend('Cluster 1','Cluster 2','Cluster 3','location','best');
title('Ground Truth Data');

a = zlim;
b = ylim;
x = max(a(2)-a(1),b(2)-b(1))/2;
xlim([-x x]);
clear a b x;
hold off;

% Print set 1
figure
hold on;
scatter(X1(cloud1,1),X1(cloud1,2),'o','blue');
scatter(X1(cloud2,1),X1(cloud2,2),'+','red');
scatter(X1(cloud3,1),X1(cloud3,2),'x','black');
legend('Cluster 1','Cluster 2','Cluster 3','location','best');
title('Dataset 1');

b = ylim;
x = (b(2)-b(1))/2;
xlim([-x x]);
clear b x;
hold off;

% Print set 2
figure
hold on;
scatter(X2(cloud1,1),X2(cloud1,2),'o','blue');
scatter(X2(cloud2,1),X2(cloud2,2),'+','red');
scatter(X2(cloud3,1),X2(cloud3,2),'x','black');
legend('Cluster 1','Cluster 2','Cluster 3','location','best');
title('Dataset 1');

b = ylim;
x = (b(2)-b(1))/2;
xlim([-x x]);
clear b x;
hold off;

% Print latent space version of set 1
figure
hold on;
myscat(P1(:,cloud1),'o','blue');
myscat(P1(:,cloud2),'+','red');
myscat(P1(:,cloud3),'x','black');
legend('Cluster 1','Cluster 2','Cluster 3','Location','best');

% Print latent space version of set 2
myscat(P2(:,cloud1),'o','blue');
myscat(P2(:,cloud2),'+','red');
myscat(P2(:,cloud3),'x','black');

b = ylim;
x = (b(2)-b(1))/2;
xlim([-x x]);
clear b x;
hold off;

% Print latent space set1 and set2 together.
figure
hold on;
scatter(P1(1,:),P1(2,:),30,[0.5 0 0.5]);
scatter(P2(1,:),P2(2,:),30,[0 0.5 0]);
legend('Set 1','Set 2','Location','best');

b = ylim;
x = (b(2)-b(1))/2;
xlim([-x x]);
clear b x;
hold off;

%% This code will plot the thing I'm calling dataset 2
% it's the one with a line of points in set1, and a line of points with an
% abrupt break in set2.

% figure
% hold on;
% scatter(X1(1:n1,1),X1(1:n1,2),'o','b');
% scatter(X1((1:n2)+n1,1),X1((1:n2)+n1,2),'+','r');
% scatter(X1((1:n3)+n1+n2,1),X1((1:n3)+n1+n2,2),'x','black');
% legend('Class 1','Class 2','Class 3','location','best');
% title('Dataset 1');
% 
% figure
% hold on;
% scatter(X2(1:n1,1),X2(1:n1,2),'o','b');
% scatter(X2((1:n2)+n1,1),X2((1:n2)+n1,2),'+','r');
% scatter(X2((1:n3)+n1+n2,1),X2((1:n3)+n1+n2,2),'x','black');
% legend('Class 1','Class 2','Class 3','location','best');
% title('Dataset 2');