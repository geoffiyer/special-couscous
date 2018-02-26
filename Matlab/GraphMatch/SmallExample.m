
theTime = zeros(1,8);
trialsPerSlot = 5;

for i = 1:size(theTime,2)
    
    [X1,X2,super] = GenerateData(1,100*2^i,100*2^i);
    % X1 = [0.0001 0.0001; 0.05 0.05; 0.1 0.1; 1 0.1; 1 0];
    % X2 = [1 1; 1.1 1.1; 2.1 2; 2 2; 2.05 2.05];
    
    super = [1 1];
    
    tic
    for j = 1:trialsPerSlot
        [U1,U2] = getEigenvectors(X1,X2,'L2');
        [signs, evecassign] = fixSigns(U1(super(:,1),:),U2(super(:,2),:));
        U2 = U2.*signs';
        U2(:, evecassign(:,2)) = U2;
        [assign1to2] = GraphMatch(U1,U2,2*round(size(U2,1)^(1/3)));
    end
    theTime(i) = toc/trialsPerSlot;
    
    %     edgeweights1to2 = dot(U1(assign1to2(:,1),:),U2(assign1to2(:,2),:),2);
    %     edgeweights2to1 = dot(U1(assign2to1(:,1),:),U2(assign2to1(:,2),:),2);
    %
    %     viewMatch(assign1to2,X1,X2,assign2to1)
    %
end

ratios = zeros(1,size(theTime,2)-1);
for i = 1:size(ratios,2)
    ratios(i) = theTime(i+1)/theTime(i);
end

% %% Look at edge weights a little
% % this code picks out the strongNum of the strongest edges, and gives you the assignment
% strongNum = 30;
% edgeweights = [assign(:,1), assign(:,2), U(sub2ind(size(U), assign(:,2), assign(:,1)))];
% edgeweights = sortrows(edgeweights,-3);
% strongAssign = edgeweights(1:strongNum,1:2);

% %% Difference between coregister and graph match
% X1Norms = sqrt(sum(X1.^2,2));
% X2Norms = sqrt(sum(X2.^2,2));
% diffNorms = zeros(size(assign1to2,1),4);
%
% diffNorms(:,1) = sqrt(sum((X1(assign1to2(:,1),:) - X1(assign1to2(:,2),:)).^2,2));
% diffNorms(:,2) = sqrt(sum((X2(assign1to2(:,1),:) - X2(assign1to2(:,2),:)).^2,2));
% diffNorms(:,3) = sqrt(sum((X1(assign2to1(:,1),:) - X1(assign2to1(:,2),:)).^2,2));
% diffNorms(:,4) = sqrt(sum((X2(assign2to1(:,1),:) - X2(assign2to1(:,2),:)).^2,2));

% % NOTE: I actually don't like the below idea. It's not invariant under
% %   isometry so it's pretty bad. But I'll save the code.
% % Calculate difference as ratio of (perm1 - perm2) / (perm1 + perm2)
% % Idea: We need some sort of normalization
% diff_vec = sqrt(2)*(X1(assign1to2(:,2),:) - X1(assign1to2(:,1),:)) ./ ...
%     (X1Norms(assign1to2(:,2)) + X1Norms(assign1to2(:,1)));
% diffNorms(:,1) = sqrt(sum(diff_vec.^2,2));
%
% diff_vec = sqrt(2)*(X2(assign1to2(:,2),:) - X2(assign1to2(:,1),:)) ./ ...
%     (X2Norms(assign1to2(:,2)) + X2Norms(assign1to2(:,1)));
% diffNorms(:,2) = sqrt(sum(diff_vec.^2,2));
%
% % % I don't understand why this line works so well
% % % We're focusing on pairs that GraphMatch says are good, but diffNorms says
% % % are bad. And somehow this picks up the most interesting data
% % diffNorms(edgeweights1to2<edgethresh1to2,1:2) = ...
% %     repmat([0 0], [size(edgeweights1to2(edgeweights1to2 < edgethresh1to2)) , 1]);
%
% diff_vec = sqrt(2)*(X1(assign2to1(:,2),:) - X1(assign2to1(:,1),:)) ./ ...
%     (X1Norms(assign2to1(:,2)) + X1Norms(assign2to1(:,1)));
% diffNorms(:,3) = sqrt(sum(diff_vec.^2,2));
%
% diff_vec = sqrt(2)*(X2(assign2to1(:,2),:) - X2(assign2to1(:,1),:)) ./ ...
%     (X2Norms(assign2to1(:,2)) + X2Norms(assign2to1(:,1)));
% diffNorms(:,4) = sqrt(sum(diff_vec.^2,2));
%
% % diffNorms(edgeweights2to1<edgethresh2to1,3:4) = ...
% %     repmat([0 0], [size(edgeweights2to1(edgeweights2to1 < edgethresh2to1)) , 1]);

% %% Plot the highlighted points
%
% xbounds = [min( [X1(:,1); X2(:,1)]) - 1, max( [X1(:,1); X2(:,1)]) + 1];
% ybounds = [min( [X1(:,2); X2(:,2)]) - 1, max( [X1(:,2); X2(:,2)]) + 1];
%
% % viewMatch(assign1to2,X1,X2,assign2to1);
% figure(1)
% subplot(1,2,1)
% scatter3(X1(:,1),X1(:,2),zeros(size(X1,1),1),'ob');
% hold on
% scatter3(X2(:,1),X2(:,2), ones(size(X2,1),1),'+r');
% line([X1(assign1to2(:,1),1) X2(assign1to2(:,2),1)]', ...
%     [X1(assign1to2(:,1),2) X2(assign1to2(:,2),2)]', [0 1],'color','blue');
% title('assign 1 to 2')
% hold off
% subplot(1,2,2)
% scatter3(X1(:,1),X1(:,2),zeros(size(X1,1),1),'ob');
% hold on
% scatter3(X2(:,1),X2(:,2), ones(size(X2,1),1),'+r');
% title('assign 2 to 1')
% line([X1(assign2to1(:,1),1) X2(assign2to1(:,2),1)]', ...
%     [X1(assign2to1(:,1),2) X2(assign2to1(:,2),2)]', [0 1],'color','red');
% hold off;
%
% % Input Data
% figure(2)
% subplot(2,2,1)
% scatter(X1(:,1),X1(:,2),[],[0 0.7 0])
% hold on
% scatter(X2(:,1),X2(:,2),[],'ob')
% xlabel('Time');
% ylabel('Observation');
% axis([xbounds ybounds]);
% legend('Modality 1','Modality 2','Location','best');
% hold off
%
% subplot(2,2,2)
% plot(diffNorms(:,1))
% hold on
% plot(diffNorms(:,2))
% title('From match 1 to 2')
% hold off
%
% subplot(2,2,3)
% plot(diffNorms(:,3))
% hold on
% plot(diffNorms(:,4))
% title('From match 2 to 1')
% hold off
%
% idx = false(size(X1,1),4);
% for i = 1:4
%     temp = FindEdgethresh(diffNorms(:,i));
%     idx(diffNorms(:,i) > temp,i) = true;
%     figure
% %     scatter(X2(idx(:,i),1),X2(idx(:,i),2),'rx');
% %     hold on
% %     scatter(X2(~idx(:,i),1),X2(~idx(:,i),2),'gx');
%     hold off
% end

% % I don't like this way of representing the data
% % Highlighted Points
% subplot(1,3,2)
% hold on
% testidx = (edgeweights1to2 <= edgethresh1to2);
% scatter(X1(testidx,1),X1(testidx,2),[],[0 0.7 0]);
% testidx = (edgeweights1to2 >= edgethresh1to2);
% scatter(X1(testidx,1),X1(testidx,2),'+r');
% testidx = (edgeweights1to2 <= edgethresh1to2);
% scatter(X2(testidx,1),X2(testidx,2),[],[0 0.7 0]);
% testidx = (edgeweights1to2 >= edgethresh1to2);
% scatter(X2(testidx,1),X2(testidx,2),'+r');
% xlabel('Time');
% ylabel('Observation');
% axis([xbounds ybounds]);
% title('Highlighted from 1 to 2 match')
% legend('boring','interesting','Location','best');
% hold off
%
% % Highlighted Points
% subplot(1,3,3)
% hold on
% testidx = (edgeweights2to1 <= edgethresh2to1);
% scatter(X1(testidx,1),X1(testidx,2),[],[0 0.7 0]);
% testidx = (edgeweights2to1 >= edgethresh2to1);
% scatter(X1(testidx,1),X1(testidx,2),'+r');
% testidx = (edgeweights2to1 <= edgethresh2to1);
% scatter(X2(testidx,1),X2(testidx,2),[],[0 0.7 0]);
% testidx = (edgeweights2to1 >= edgethresh2to1);
% scatter(X2(testidx,1),X2(testidx,2),'+r');
% xlabel('Time');
% ylabel('Observation');
% axis([xbounds ybounds]);
% title('Highlighted from 2 to 1 match')
% legend('boring','interesting','Location','best');
% hold off

% figure
% hold on
% idx = (edgeweights1to2 > edgethresh1to2);
% scatter(X1(idx,1),X1(idx,2),'ob');
% scatter(X1(~idx,1),X1(~idx,2),'black','x');
% idx = (edgeweights1to2 > edgethresh1to2);
% scatter(X2(idx,1),X2(idx,2),'+r');
% scatter(X2(~idx,1),X2(~idx,2),'black','x');
% % line([X1(assign1to2(:,1),1) X2(assign1to2(:,2),1)]', ...
% %        [X1(assign1to2(:,1),2) X2(assign1to2(:,2),2)]', [0 1],'color','blue');
%
% figure
% hold on
% idx = (edgeweights2to1 > edgethresh2to1);
% scatter3(X1(idx,1),X1(idx,2),zeros(size(X1(idx,1),1),1),'ob');
% scatter3(X1(~idx,1),X1(~idx,2),zeros(size(X1(~idx,1),1),1),'black','x');
% idx = (edgeweights2to1 > edgethresh2to1);
% scatter3(X2(idx,1),X2(idx,2), ones(size(X2(idx,1),1),1),'+r');
% scatter3(X2(~idx,1),X2(~idx,2), ones(size(X2(~idx,1),1),1),'black','x');
% line([X1(assign2to1(:,1),1) X2(assign2to1(:,2),1)]', ...
%         [X1(assign2to1(:,1),2) X2(assign2to1(:,2),2)]', [0 1],'color','blue');