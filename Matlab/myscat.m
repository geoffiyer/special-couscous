function [] = myscat( X , label, color)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if nargin < 3
    label = 'o';
    color = 'blue';
end

if( size(X,1) == 2)
    scatter(X(1,:),X(2,:), 80 ,label, color);
elseif (size(X,1) == 3)
    scatter3(X(1,:),X(2,:),X(3,:), 80, label, color);
end

end

