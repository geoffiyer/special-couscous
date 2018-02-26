
n1 = 20;    % size of main cloud
n2 = 5;     % size of secondary cloud
xmin = 0;
xmax = 20;
ysdev = 1;
ymean = 0;

cloud2ymean = 7;
cloud2xmin = 15;
cloud2xmax = 18;

X1 = zeros(2,n1+n2);
X2 = zeros(2,n1+n2);
X1(1,1:n1) = rand(1,n1)*(xmax-xmin) + xmin;
X1(1,(1+n1):(n1+n2)) = rand(1,n2)*(cloud2xmax - cloud2xmin) + cloud2xmin;
X2(1,:) = X1(1,:);

X1(2,:) = randn(1,n1+n2)*ysdev + ymean;
X2(2,:) = X1(2,:);
X2(2,(1+n1):(n1+n2)) = X2(2,(1+n1):(n1+n2)) + (cloud2ymean - ymean);