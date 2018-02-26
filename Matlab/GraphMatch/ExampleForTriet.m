%% Generate Data and do graph match

% I have the hungarian algorithm match working pretty well,
% and the naive max match sort of working.
% True means use hungarian, false means use naive max
useHungarian = false;

[X1,X2,super,perm] = GenerateData(3, 1000, 1000, false);
X2(:,2) = X2(:,2) + 10;

[U1,U2] = getEigenvectors(X1,X2,'L2');
[signs, evecassign] = fixSigns(U1(super(:,1),:),U2(super(:,2),:));
U2 = U2.*signs';
U2(:, evecassign(:,2)) = U2;

if(useHungarian)
    [assign,Y1,Y2,Yassign,U1,U2] = GraphMatch(X1,X2,super);
    assign = sortrows(assign,1);
else
    [assign1to2,assign2to1] = GraphMatch2(U1,U2,1);
    assign = sortrows(assign1to2,1);
end
% % Leftover from when we still calculated the full matrix U
% edgeweights = U(sub2ind(size(U), assign(:,1), assign(:,2)));
edgeweights = dot(U1(assign(:,1),:),U2(assign(:,2),:),2);
edgethresh = FindEdgethresh(edgeweights);

%% Difference between coregister and graph match
X1Norms = sqrt(sum(X1.^2,2));
X2Norms = sqrt(sum(X2.^2,2));
diffNorms = zeros(size(assign,1),4);

% Calculate difference as ratio of (perm1 - perm2) / (perm1 + perm2)
% Idea: We need some sort of normalization
diff_vec = 2*(X1(assign(:,2),:) - X1(assign(:,1),:)) ./ (X1Norms(assign(:,2)) + X1Norms(assign(:,1)));
diffNorms(:,1) = sqrt(sum(diff_vec.^2,2));

diff_vec = 2*(X2(assign(:,2),:) - X2(assign(:,1),:)) ./ (X2Norms(assign(:,2)) + X2Norms(assign(:,1)));
diffNorms(:,2) = sqrt(sum(diff_vec.^2,2));

% % I don't understand why this line works so well
% % We're focusing on pairs that GraphMatch says are good, but diffNorms says
% % are bad. And somehow this picks up the most interesting data
% diffNorms(edgeweights<edgethresh,1:2) = repmat([0 0], [size(edgeweights(edgeweights < edgethresh)) , 1]);

if(useHungarian)
    assign = sortrows(assign,2);
else
    assign = sortrows(assign2to1,2);
end

diff_vec = sqrt(2)*(X1(assign(:,2),:) - X1(assign(:,1),:)) ./ (X1Norms(assign(:,2)) + X1Norms(assign(:,1)));
diffNorms(:,3) = sqrt(sum(diff_vec.^2,2));

diff_vec = sqrt(2)*(X2(assign(:,2),:) - X2(assign(:,1),:)) ./ (X2Norms(assign(:,2)) + X2Norms(assign(:,1)));
diffNorms(:,4) = sqrt(sum(diff_vec.^2,2));

% diffNorms(edgeweights<edgethresh,3:4) = repmat([0 0], [size(edgeweights(edgeweights < edgethresh)) , 1]);

%% Plot some stuff

xbounds = [0,20];
ybounds = [min( [X1(:,2); X2(:,2)]) - 1, max( [X1(:,2); X2(:,2)]) + 1];

% Input Data
fig = figure;
hold on
scatter(X1(:,1),X1(:,2),50,[0 0.7 0])
scatter(X2(:,1),X2(:,2),50,'ob')
xlabel('Time');
ylabel('Observation');
axis([xbounds ybounds]);
legend('Modality 1','Modality 2','Location','best');
set(gca,'fontsize',25)
saveas(fig,'/home/gsiyer/Schoolwork/Chanussot/Presentations/GISPA_2017-10-19/Images/ChangeDetect/data.png');
hold off

% Points by distance thing
fig = figure;
hold on
diffnorms = sqrt(sum((X1(assign2to1(:,1),:) - X1(assign2to1(:,2),:)).^2,2));
diffnorms = diffnorms/max(max(diffnorms));
scatter(X1(:,1),X1(:,2),40,diffnorms,'filled');
diffnorms = sqrt(sum((X2(assign2to1(:,1),:) - X2(assign2to1(:,2),:)).^2,2));
diffnorms = diffnorms/max(max(diffnorms));
scatter(X2(:,1),X2(:,2),40,diffnorms,'filled');
colorbar;
xlabel('Time');
ylabel('Observation');
axis([xbounds ybounds]);
set(gca,'fontsize',25)
hold off
saveas(fig,'/home/gsiyer/Schoolwork/Chanussot/Presentations/GISPA_2017-10-19/Images/ChangeDetect/norms.png');

% % Highlighted Points
% subplot(1,2,2)
% hold on
% testidx = (edgeweights < edgethresh);
% scatter(X1(testidx,1),X1(testidx,2),[],[0,0.7,0]);
% testidx = (edgeweights < edgethresh);
% scatter(X2(testidx,1),X2(testidx,2),[],'ob');
% testidx = (edgeweights >= edgethresh);
% scatter(X1(testidx,1),X1(testidx,2),'+r');
% testidx = (edgeweights >= edgethresh );
% scatter(X2(testidx,1),X2(testidx,2),'+r');
% xlabel('Time');
% ylabel('Observation');
% axis([xbounds ybounds]);
% legend('Modality 1','Modality 2','Highlighted Matches','Location','best');
% hold off
