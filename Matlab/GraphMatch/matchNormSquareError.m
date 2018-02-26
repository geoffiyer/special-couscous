function [ error ] = matchNormSquareError( assign, W1, W2 )
%MATCHERROR Summary of this function goes here
%   Detailed explanation goes here

assign = sortrows(assign,1);
newW2 = W2(:,assign(:,2));
newW2 = newW2(assign(:,2),:);

error = norm(W1 - newW2,'fro')^2;

end

