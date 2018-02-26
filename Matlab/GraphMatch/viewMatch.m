function [] = viewMatch(assign1to2,X1,X2,assign2to1)
%% Here X1 and X2 are datasets
% assign is a n x 2 array, where n is number of matches.
% assign(i,1) and assign(i,2) are the matched pair in X1, X2 (respectively)

% Lately I've been writing code where there aren't always one-to-one
% matches. So I've introduces assignType and assign2to1. Here's how it
% works: if assignType == 'bijection', then do the 1-to-1 thing. If
% assignType == 'anything else here', then just break the whole thing
% because I don't want to write that code yet.


if(~exist('assign2to1','var'))
    figure
    hold on
    scatter3(X1(assign1to2(:,1),1),X1(assign1to2(:,1),2),zeros(size(X1,1),1),'ob');
    scatter3(X2(assign1to2(:,2),1),X2(assign1to2(:,2),2),ones(size(X2,1),1),'+r');
    line([X1(assign1to2(:,1),1) X2(assign1to2(:,2),1)]', [X1(assign1to2(:,1),2) X2(assign1to2(:,2),2)]', [0 1]);
    hold off;
else
    figure
    subplot(1,2,1)
    hold on
    scatter3(X1(:,1),X1(:,2),zeros(size(X1,1),1),'ob');
    scatter3(X2(:,1),X2(:,2), ones(size(X2,1),1),'+r');
    line([X1(assign1to2(:,1),1) X2(assign1to2(:,2),1)]', ...
        [X1(assign1to2(:,1),2) X2(assign1to2(:,2),2)]', [0 1],'color','blue');
    title('assign 1 to 2')
    subplot(1,2,2)
    hold on
    scatter3(X1(:,1),X1(:,2),zeros(size(X1,1),1),'ob');
    scatter3(X2(:,1),X2(:,2), ones(size(X2,1),1),'+r');
    title('assign 2 to 1')
    line([X1(assign2to1(:,1),1) X2(assign2to1(:,2),1)]', ...
        [X1(assign2to1(:,1),2) X2(assign2to1(:,2),2)]', [0 1],'color','red');
    hold off;
end

end

% %% This was some code about viewing 3-d to 1-d matches.
% % I feel like I should delete it forever but instead I'll leave it here
% figure
% hold on
% scatter3(X1(assign(:,1),1),X1(assign(:,1),2),X1(assign(:,1),3),'ob');
% scatter3(X2(assign(:,2),1),zeros(size(assign,1),1),zeros(size(assign,1),1),'+r'); 
% line([X1(assign(:,1),1) X2(assign(:,2),1)]', [X1(assign(:,1),2) zeros(size(assign,1),1)]', [X1(assign(:,1),3) zeros(size(assign,1),1)]);


