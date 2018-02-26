% Plan: Take data from the GenerateData function. Assume that I don't know
% any class labels, but I do know the point correspondences between X_1 and
% X_2. Try to reconstruct orignal figure using the similarity matrix idea.

% Creates 3 clouds of 3 points each
% cloud_1 centered at [0,0,0], it's the red +
% cloud_2 centered at [0,20,0], it's the blue o
% cloud_3 centered at [0,0,20], it's the green x
% in X_1, we smash cloud_3 into cloud_1
% in X_2, we smash cloud_2 into cloud_1
GenerateData;
X_Full = [X_1, X_2];
cloud_1 = X_Full(3,:)==1;
cloud_2 = X_Full(3,:)==2;
cloud_3 = X_Full(3,:)==3;

W_size = size(X_Full,2);
W_s = zeros(W_size,W_size);
for i = 1:W_size
    for j = (i+1):W_size
        if( j == i+n*num_clouds)
            W_s(i,j) = 1;
            W_s(j,i) = 1;
        end
    end
end
L_s = diag(sum(W_s,2)) - W_s;

W = zeros(W_size);
for i = 1:W_size
    if(i <= n*num_clouds)
        for j = 1:n*num_clouds
            W(i,j) = exp(-1*norm(X_Full(1:2,i) - X_Full(1:2,j))^2);
        end
    else
        for j = (n*num_clouds+1):2*n*num_clouds
            W(i,j) = exp(-1*norm(X_Full(1:2,i) - X_Full(1:2,j))^2);
        end
    end
end
L = diag(sum(W,2)) - W;

[V,D] = eig(L+L_s);
I = find(diag(D)>100*eps,2);
if( I(1) ~= 2 )
    warning('Weird eigenvalue stuff. Check D');
end

figure;
hold on;
scatter(V(cloud_1,I(1)),V(cloud_1,I(2)),'+','red');
scatter(V(cloud_2,I(1)),V(cloud_2,I(2)),'o','blue');
scatter(V(cloud_3,I(1)),V(cloud_3,I(2)),'x','green');
hold off;

figure;
hold on;
scatter3(X_G(1,cloud_1(1:9)),X_G(2,cloud_1(1:9)),X_G(3,cloud_1(1:9)),'+','red');
scatter3(X_G(1,cloud_2(1:9)),X_G(2,cloud_2(1:9)),X_G(3,cloud_2(1:9)),'o','blue');
scatter3(X_G(1,cloud_3(1:9)),X_G(2,cloud_3(1:9)),X_G(3,cloud_3(1:9)),'x','green');
hold off;


% figure
% hold on;
% scatter(X_1(1,cloud_1(1:9)),X_1(2,cloud_1(1:9)),'+','red');
% scatter(X_1(1,cloud_2(1:9)),X_1(2,cloud_2(1:9)),'o','blue');
% scatter(X_1(1,cloud_3(1:9)),X_1(2,cloud_3(1:9)),'x','green');
% xlim([-10 10])
% hold off;
% 
% figure
% hold on;
% scatter(X_2(1,cloud_1(1:9)),X_2(2,cloud_1(1:9)),'+','red');
% scatter(X_2(1,cloud_2(1:9)),X_2(2,cloud_2(1:9)),'o','blue');
% scatter(X_2(1,cloud_3(1:9)),X_2(2,cloud_3(1:9)),'x','green');
% xlim([-10 10])
% hold off;