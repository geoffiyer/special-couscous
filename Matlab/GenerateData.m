
%% 3 clouds in R^3, projected in different ways to make X1 and X2

n = 100;  % Number of data points per cloud.
num_clouds = 3;
% Centers of our 3 gaussian clouds
mu = [0, 0, 0 ; 0,20,0 ; 0,0,20];

% Ground truth. Data in R^3
% Each individual dimension is random normal.
% Do this instead of multivariate random normal
% because apparently that costs $200 for the package
XG = randn(3,num_clouds*n) + repmat(mu,1,n);
XG = [XG ; repmat(1:num_clouds,1,n)];  %4th dim stores the class

% This code hasn't been updated for num_clouds
X1 = XG([1,2,4],:);   % Project onto the xy plane
                        % Cloud 3 merges with cloud 1
X2 = XG([1,3,4],:);   % Project onto the xz place
                        % Cloud 2 merges with cloud 1
                        
cloud1 = 1:3:(3*n);
cloud2 = 2:3:(3*n);
cloud3 = 3:3:(3*n);
