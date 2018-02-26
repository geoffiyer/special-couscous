
GenerateData;
X1 = X1/6;
X2 = X2/6;

[nrx1, ncx1] = size(X1);
[nrx2, ncx2] = size(X2);

W = zeros(ncx1,ncx1);

for i = 1:ncx1
    for j = (i+1):ncx1
        W(i,j) = exp(-1*max( norm(X1(:,i)-X1(:,j)), norm(X2(:,i)-X2(:,j))));
        W(j,i) = W(i,j);
    end
end

L = diag(sum(W,2)) - W;

[V,D] = eig(L);

figure
hold on;
myscat(V(cloud1,2:3)','+','red');
myscat(V(cloud2,2:3)','o','blue');
myscat(V(cloud3,2:3)','x','black');
h_legend = legend('Cluster 1','Cluster 2','Cluster 3','Location','best');
set(h_legend,'FontSize',26);
xlabel('Eigenvector 1') % x-axis label
ylabel('Eigenvector 2') % y-axis label
hold off;

