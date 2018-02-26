function [ d ] = angledist(a,b )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

d = 1-((a*b')./(sqrt(sum(a.^2,2))*sqrt(sum(b.^2,2))'));
d = sqrt(d);

%% Slow way with for loop that I used for FAR too long.
% d = zeros(size(a,1),size(b,1));
% 
% for i = 1:size(a,1)
%     num1 = norm(a(i,:));
%     for j = 1:size(b,1)
%         d(i,j) = 1-dot(a(i,:),b(j,:)) / (num1 * norm(b(j,:)));
%     end
% end
% 
% d = sqrt(d);

end