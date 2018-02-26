function [ E ] = permEnergy( p, X_1, X_2 )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
E = 0;

% sqrt(sum((repmat(X_1(1:2,i),1,3)-X_1(1:2,1:3)).^2,1))
% gives distance from chosen point i to each preset correspondence point j

for i = 4:9
    E = E + sum( max( ...
        sqrt(sum((repmat(X_1(1:2,i),1,3)-X_1(1:2,1:3)).^2,1)), ...
        sqrt(sum((repmat(X_2(1:2,p(i-3)),1,3)-X_2(1:2,1:3)).^2,1)) ...
        ) );
end

end

