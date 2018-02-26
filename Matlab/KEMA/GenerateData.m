
if set1
    
    n = 80;  % Number of data points per cloud.
    num_clouds = 3;
    % Centers of our 3 gaussian clouds
    mu = [0, 0, 0 ; 0,20,0 ; 0,0,20];
    
    % Ground truth. Data in R^3
    % Each individual dimension is random normal.
    % Do this instead of multivariate random normal
    % because apparently that costs $200 for the package
    XG = randn(num_clouds*n,3) + repmat(mu,n,1);
    XG = [XG  repmat((1:num_clouds)',n,1)];  %4th dim stores the class
    
    % This code hasn't been updated for num_clouds
    X1 = XG(:,[1,2]);   % Project onto the xy plane
    % Cloud 3 merges with cloud 1
    X2 = XG(:,[1,3]);   % Project onto the xz place
    % Cloud 2 merges with cloud 1
    
    X = [X1 ; X2];
    Y = [zeros(size(X1,1),1)+1 ; zeros(size(X2,1),1)+2];
    Y1 = XG(:,4);
    Y2 = XG(:,4);
    
    cloud1 = 1:3:(3*n);
    cloud2 = 2:3:(3*n);
    cloud3 = 3:3:(3*n);
    
    N = 20;
    r = 1;
    
elseif set2
    
    n1 = 200;    % size of main cloud
    n2 = 200;     % size of secondary cloud
    n3 = 200;
    n = n1+n2+n3;
    xmin = 0;
    xmax = 20;
    ysdev = 1;
    ymean = 0;
    
    cloud2ymean = 7;
    cloud2xmin = 15;
    cloud2xmax = 18;
    
    X1 = zeros(n,2);
    X1(1:n,1) = rand(n,1)*(xmax-xmin) + xmin;
    X1 = sort(X1);
    X1(:,2) = randn(n,1)*ysdev + ymean;
    X2 = X1;
    X2((1+n1+n2):n,2) = X2((1+n1+n2):n,2) + (cloud2ymean - ymean);
    
    X = [X1 ; X2];
    Y = [ones(n,1);zeros(n,1)+2];
    Y1 = [ones(n1,1);ones(n2,1)+1;ones(n3,1)+2];
    Y2 = Y1;
    
    clear cloud2xmax cloud2xmin cloud2ymean n n1 n2 n3 xmax xmin ymean ysdev
    
elseif set3
    
    n1 = 200;
    n2 = 200;
    n3 = 200;
    n = n1+n2+n3;
    line1 = [20 0];
    cloud1 = [4 -8];
    sdev = 0.5;
    line2 = [-10 20;];
    cloud2 = [8 0];
        
    X1 = zeros(n,2);
    X1((1+n1):n,:) = rand(n3+n2,1) * line1 + randn(n3+n2,2)*sdev;
    X1((1+n1):n,:) = sortrows(X1((1+n1):n,:));
    X1(1:n1,:) = randn(n1,2)*sdev + cloud1;
    
    X2 = zeros(n,2);
    temp = sort(rand(n3+n2,1))*(pi/2);
    X2((1+n1):n,1) = sin(temp)*-8+8;
    X2((1+n1):n,2) = cos(temp)*-8;
    X2((1+n1):n,:) = X2((1+n1):n,:) + randn(n3+n2,2)*sdev;
    X2(1:n1,:) = randn(n1,2)*sdev + cloud2;
    
    figure
    scatter(X1(:,1),X1(:,2),30,'ob');
    figure
    scatter(X2(:,1),X2(:,2),30,'+r');
    
    X = [X1 ; X2];
    Y = [ones(n,1);zeros(n,1)+2];
    Y1 = [ones(n1,1);ones(n2,1)+1;ones(n3,1)+2];
    Y2 = Y1;
    
elseif set4
    
    n = 2000;
    tmin = -pi/2;
    tmax = 2*pi + pi/2;
    tvals = (tmin:((tmax - tmin)/(n-1)):tmax)';
    
    X1 = [(tvals - tmin + 2), (tvals - tmin + 2), ones(n,1)] .* [cos(tvals),sin(tvals),tvals] + ([1, 1, 0.2] .*  rand(n,3)*2) ;
    
    % clear tvals tmin tmax n;

    scatter(X1(:,1),X1(:,2),20,'+r')
    
end
