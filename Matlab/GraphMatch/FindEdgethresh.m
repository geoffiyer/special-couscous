function [thresh, center] = FindEdgethresh(weights)

addpath ./OldNystrom

[~,center,~] = eff_kmeans(weights,2,100);
thresh = (center(1) + center(2))/2;
center = sort(center);

rmpath ./OldNystrom

% % This seems like a nice idea: find threshold based on maximum contrast.
% % But in practice I get a bunch of weirdo answers
%     sort(weights);
%     n = size(weights,1);
%     delta = round(max(n/500,5));
%     
%     maxContrast = 0;
%     idx = 1;
%     
%     for i=1:(n-delta)
%         if(weights(min(i+delta,n)) - weights(i) > maxContrast*1.1 && i > n/3)
%             maxContrast = weights(i+delta) - weights(i);
%             idx = i;
%         end
%     end
%     
%     thresh = weights(idx);

end