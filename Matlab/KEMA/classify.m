% class = classify(sample,training,group) 
% classifies each row of the data in sample into one of the groups in
% training. sample and training must be matrices with the same number of
% columns. group is a grouping variable for training. Its unique values
% define groups; each element defines the group to which the corresponding
% row oftraining belongs. group can be a categorical variable, a numeric
% vector, a character array, or a cell array of character vectors. training
% and group must have the same number of rows. classify treats NaNs or
% empty character vectors in group as missing values, and ignores the
% corresponding rows of training. The output class indicates the group to
% which each row of sample has been assigned, and is of the same type as
% group.

function [class,P] = classify(sample, training, group)

W = LDA(training, group);
L = [ones(size(sample,1),1) sample] * W';
P = exp(L) ./ repmat(sum(exp(L),2),[1 size(W,1)]);
[~,class] = max(P,[],2);

end