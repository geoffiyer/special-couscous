
GenerateData;

X_Full = [X_1, X_2];
cloud_1 = X_Full(3,:)==1;
cloud_2 = X_Full(3,:)==2;
cloud_3 = X_Full(3,:)==3;
% Remove some labels!
X_Full(3,[1,5,9,13,17]) = [0,0,0,0,0];

W_size = size(X_Full,2);
W_s = zeros(W_size,W_size);
for i = 1:W_size
    for j = 1:W_size
        if( X_Full(3,i) == X_Full(3,j) && X_Full(3,i) ~= 0 )
            W_s(i,j) = 2;
        end
        %if( X_Full(3,i) ~= X_Full(3,j) && X_Full(3,i) ~= 0 && X_Full(3,j) ~= 0)
        %    W_s(i,j) = 0;
        %end
        % W(i,j) = W(i,j) + norm(X_Full(1:2,i)-X_Full(1:2,j));
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
find(D>0.001,2);

h = figure;
hold on;
scatter(V(cloud_1,2),V(cloud_1,3),'+','red');
scatter(V(cloud_2,2),V(cloud_2,3),'o','blue');
scatter(V(cloud_3,2),V(cloud_3,3),'x','green');
hold off;

% [V,D] = eig(L);
% Latent = V(:,(W_size-1):W_size)';
% addpath ./kmeans/
% Labels = kmeans(Latent,3);
% rmpath ./kmeans/