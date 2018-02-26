function [ K1 ] = matchClasses( K1,K2 )
% Input two classifications. Renumber classes in one to match the other. 
%   INPUT:
%      2d arrays  K1, K2 representing classifications of the same image.
%        expect size(K1) == size(K2), and unique(K1) == unique(K2)
%   Try to permute the values of K1 so that the numbering of classes
%   matches with K2.
%   Output:
%      Edited K1

% Actual algorithm: Greedy based on biggest class.

assert(isequal(size(K1),size(K2)), 'K1 and K2 must be same size');
numClasses = max(max(max(K1)),max(max(K2)));

classSizes = zeros(numClasses,1);
for i = 1:size(classSizes,1)
    classSizes(i) = sum(sum(K1 == i));
end
[~,idx] = sort(classSizes,'descend');

perm1 = zeros(numClasses,1);
perm2 = (1:numClasses)';

for i = 1:size(perm1,1)
    bestMatch = 0;
    bestIdx = 0;
    for j = 1:size(perm2,1)
        temp = sum(sum( (K1 == idx(i)) & (K2 == perm2(j)) ));
        if(temp >= bestMatch)
            bestMatch = temp;
            bestIdx = j;
        end
    end
    perm1(i) = perm2(bestIdx);
    perm2(bestIdx) = [];
end

tempK1 = K1;
for i = 1:size(perm1,1)
    K1(tempK1 == idx(i)) = perm1(i);
end

