
%% Compare coarse graph to full graph (debugging some indexing stuff)
% 
% for i = 1:size(Y2,1)
% figure
% hold on
% scatter(X2(:,1),X2(:,2),'xg');
% scatter(X2(clusters2(i,:),1),X2(clusters2(i,:),2),'ob');
% scatter(Y2(i,1),Y2(i,2),'+r');
% hold off
% end

i = 1;
assign = zeros(clusterSize,2);
idx1 = clusters1(smallAssign(i,1),:);
idx2 = clusters2(smallAssign(i,2),:);
Usmall = U(idx2,idx1);
[inverse_assign,~] = munkres(1 - Usmall);
local_assign = [inverse_assign', (1:size(inverse_assign,2))'];
assign(((i-1)*clusterSize+1):i*clusterSize,:) = [clusters1(smallAssign(i,1),local_assign(:,1))', ...
    clusters2(smallAssign(i,2),local_assign(:,2))'];

figure
hold on;
scatter3(X1(:,1),X1(:,2),zeros(size(X1,1),1),'xg');
scatter3(X2(:,1),X2(:,2),ones(size(X2,1),1),'xg');
scatter3(X1(clusters1(smallAssign(i,1),:),1),X1(clusters1(smallAssign(i,1),:),2),zeros(clusterSize,1),'ob');
scatter3(X2(clusters2(smallAssign(i,2),:),1),X2(clusters2(smallAssign(i,2),:),2),ones(clusterSize,1),'+r');
line([X1(assign(:,1),1) X2(assign(:,2),1)]', [X1(assign(:,1),2) X2(assign(:,2),2)]', [0 1]);
hold off;

%% Plotting some stuff for Jocelyn

% Assumes X1,Y1,X2,Y2 have already been defined

% figure
% subplot(1,2,1)
% scatter(X2(:,1),X2(:,2))
% title('Original Data')
% subplot(1,2,2)
% scatter(Y2(:,1),Y2(:,2))
% title('Coarsified Data')


figure
subplot(2,3,1)
scatter(X1(:,1),X1(:,2),'ob')
title('Data 1')
subplot(2,3,2)
scatter(X2(:,1),X2(:,2),'+r')
title('Data 2')
subplot(2,3,3)
scatter3(X1(:,1),X1(:,2),zeros(size(X1,1),1),'ob');
hold on;
scatter3(X2(:,1),X2(:,2),ones(size(X2,1),1),'+r');
line([X1(assign(:,1),1) X2(assign(:,2),1)]', [X1(assign(:,1),2) X2(assign(:,2),2)]', [0 1]);
title('Full graph matching (unsupervised)')
hold off;
subplot(2,3,4)
scatter(Y1(:,1),Y1(:,2),'ob')
title('Coarsified Data 1')
subplot(2,3,5)
scatter(Y2(:,1),Y2(:,2),'+r')
title('Coarsified Data 2')
subplot(2,3,6)
scatter3(Y1(:,1),Y1(:,2),zeros(size(Y1,1),1),'ob');
hold on
scatter3(Y2(:,1),Y2(:,2),ones(size(Y2,1),1),'+r');
line([Y1(smallAssign(:,1),1) Y2(smallAssign(:,2),1)]', [Y1(smallAssign(:,1),2) Y2(smallAssign(:,2),2)]', [0 1]);
title('Coarse graph matching (unsupervised)')
hold off;