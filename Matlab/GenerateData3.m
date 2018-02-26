
    n1 = 200;    % size of main cloud
    n2 = 100;     % size of secondary cloud
    n3 = 300;
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
    
    m = 1+n1+n2;
    l = min(X2(m:n,1));
    X2(m:2:n,2) = X2(m:2:n,2) + (X2(m:2:n,1)-l)*2;

    figure
    scatter(X1(:,1),X1(:,2));
    ylim([-5 20])
    figure
    scatter(X2(:,1),X2(:,2),'+r');
    
    clear n1 n2 n3 n xmin xmax ysdev ymean cloud2ymean cloud2xmin cloud2xmax m l;
%     X = [X1 ; X2];
%     Y = [ones(n,1);zeros(n,1)+2];
%     Y1 = [ones(n1,1);ones(n2,1)+1;ones(n3,1)+2];
%     Y2 = Y1;