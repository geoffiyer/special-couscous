function [newK] = matchKtoGT( K,GT )
% Input two classifications. Renumber classes in one to match the other. 
%   INPUT:
%      2d arrays  K, representing classifications of the same image.
%                GT, representing classifications of the same image.
%        expect size(K) == size(GT). unique(GT) almost equals unique(K),
%        except that GT also has the 0 class to mean unclassified
%
%   Try to permute the values of K so that the numbering of classes
%      matches with GT.
%   Output:
%      Edited K

% Actual algorithm: Hungarian algorithm to maximize total number of
% agreeing pixels

assert(isequal(size(K),size(GT)), 'K and GT must be same size');

numClasses = max(max(K(:)),max(GT(:)));

agreemat = zeros(numClasses);

for i = 1:numClasses
    for j = 1:numClasses
        agreemat(i,j) = sum(sum(K == i & GT == j));
    end
end

addpath /home/gsiyer/Schoolwork/Chanussot/Matlab/GraphMatch/

assign = Hungarian(agreemat);

newK = zeros(size(K));

for i = 1:numClasses
    newK(K == i) = assign(i,2);
end