
function [X1,X2,super,perm] = GenerateData(version, n, m, permute)
%% Generate data from 5 different methods
% version = 1, 2, 3, 4, or 5
% n = number of points in set 1
% m = number of points in set 2
% permute = should I permute set 2 (for testing purposes)?

% Note about permutation:
% X1(perm(1),:) matches X2(1,:)
% super is edited to handle permutation
%   i.e. X1(super(1,1),:) matches X2(super(1,2),:)
% Also, I actually haven't tested super and perm since the last edit
% so who knows really...

if(nargin < 4)
    permute = false;
end

if(nargin < 3)
    m = 100;
end

if(nargin < 2)
    n = 100;
end

if( version == 1)
    
    scaling = 40;
    
    X1 = randn(2*n,2) + scaling * [repmat([1 0],n,1); repmat([-1 0],n,1)];
    X2 = randn(2*m,2) + scaling * [repmat([1 0],m,1); repmat([-1 0],m,1)];
    
    super = [(1:10:2*n)',(1:(10*m/n):2*m)'];
    
elseif(version == 2)
    
    line1 = [20 0];
    sdev = 0.5;
    
    X1 = zeros(n,2);
    X1(1:n,:) = rand(n,1) * line1 + randn(n,2)*sdev;
    X1(1:n,:) = sortrows(X1);
    
    X2 = zeros(m,2);
    temp = sort(rand(m,1))*(pi/2);
    X2(1:m,1) = cos(temp)*-8;
    X2(1:m,2) = sin(temp)*-8+8;
    X2(1:m,:) = X2(1:m,:) + randn(m,2)*sdev;
    
    % temp = round(1:(n1/m1):n1);
    % full_super = [temp' , (1:size(temp,2))'];
    % temp = round((n1+1):((n2+n3)/(m2+m3)):n);
    % full_super = [ full_super ; temp' , ((m1+1):(m1+size(temp,2)))'];
    
    full_super = [ceil(1:n/m:n)'  (1:m)'];
    super = full_super(1:10:size(full_super,1),:);
    super = [super; n m];
    
elseif(version == 3)
    
    xmin = 0;
    xmax = 20;
    ysdev = 1;
    ymean = 0;
    
    X1 = zeros(n,2);
    X1(1:n,1) = rand(n,1)*(xmax-xmin) + xmin;
    X1 = sort(X1);
    X1(:,2) = randn(n,1)*ysdev + ymean;
    
    cutpoint = round(m/2);
    
    X2 = zeros(n,2);
    X2(1:m,1) = rand(m,1)*(xmax-xmin) + xmin;
    X2 = sort(X2);
    X2(:,2) = randn(m,1)*ysdev + ymean;
    l = min(X2((1+cutpoint):m,1));
    X2((1+cutpoint):2:m,2) = X2((1+cutpoint):2:m,2) + (X2((1+cutpoint):2:m,1)-l)/(1-cutpoint/m);
    
    super = [(2:(10*n/m):n)',(2:10:m)'];
    
elseif(version == 4)
    
    n = round(n/3)*3;
    m = round(m/3)*3;
    n1 = n/3; m1 = m/3;
    n2 = n/3; m2 = m/3;
    n3 = n/3; m3 = m/3;
    line1 = [20 0];
    cloud1 = [4 -8];
    sdev = 0.5;
    cloud2 = [8 0];
    
    X1 = zeros(n,2);
    X1((1+n1):n,:) = rand(n3+n2,1) * line1 + randn(n3+n2,2)*sdev;
    X1((1+n1):n,:) = sortrows(X1((1+n1):n,:));
    X1(1:n1,:) = randn(n1,2)*sdev + cloud1;
    
    X2 = zeros(m,2);
    temp = sort(rand(m3+m2,1))*(pi/2);
    X2((1+m1):m,1) = sin(temp)*-8+8;
    X2((1+m1):m,2) = cos(temp)*-8;
    X2((1+m1):m,:) = X2((1+m1):m,:) + randn(m3+m2,2)*sdev;
    X2(1:m1,:) = randn(m1,2)*sdev + cloud2;
    
    super = [(1:10:m1)',(1:10:m1)';(n1+1:10:n1+m2)',(m1+1:10:m1+m2)';(n1+n2+1:10:n1+n2+m1)',(m1+m2+1:10:m)';];
    super = [super; n , m];
    
elseif(version == 5)
    
    X2 = [rand(m/3,2); rand(m/3,2) + repmat([10,0],m/3,1)];
    X1 = [rand(n/3,2); rand(n/3,2) + repmat([10,0],n/3,1);  ...
        rand(n/3,2) + repmat([5,8],n/3,1)];
    
    super = [ (1:(10*n/m):(2*n/3))' , (1:10:m)'];
    
elseif(version == 6)
    
    X1 = zeros(n,m,3);
    X2 = zeros(n,m,3);
    center = round([(n+1)/2,(m+1)/2]);
    
    changeMat = [ 0.8  0.8  0.4; ...
                  0.4  0.8  0.8; ...
                  0.8  0.4  0.8];
    
    for i = 1:n
        for j = 1:m
            X1(i,j,:) = (i+j)/(n+m) * [0 1 0] + (n-i+m-j)/(n+m) * [1 0 0];
            
            if( norm( [i,j] - center) < min(n,m)/6 )
                X2(i,j,:) = changeMat*[0; 0; 1];
            else
                temp = reshape(X1(i,j,:),[3 1]);
                X2(i,j,:) = changeMat*temp;
            end
        end
    end
    
    super = [ (1:(n*m))', (1:(n*m))'];

elseif (version == 7)
    
    ysdev = 0.2;
    ymean = 0;
    
    X1 = zeros(n,2);
    X1(1:n,1) = rand(n,1)*ysdev;
    X1(:,2) = randn(n,1)*ysdev + ymean;
    
    cutpoint = round(3*m/4);
    
    X2 = zeros(m,2);
    X2(1:m,1) = rand(m,1)*ysdev;
    X2(:,2) = randn(m,1)*ysdev + ymean;
    X2(cutpoint:m,2) = X2(cutpoint:m,2) + 10;
    X2(:,2) = X2(:,2) + 5;
    
    super = [(2:(10*n/m):n)',(2:10:m)'];
    
end

if(permute)
    
    % permute X2, because that's what we'd expect in the real world
    perm = randperm(size(X2,1));
    invperm = zeros(1,size(X2,1));
    invperm(perm) = 1:size(X2,1);
    X2 = X2(perm,:);
    super(:,2) = invperm(super(:,2))';
    
else
    perm = 1:size(X2);
end

end