function [ labels] = connComp( K )
% K is a 2D matrix of class values (integers). Find connected components

% Note: I use 8-way adjacency.

labels = zeros(size(K));
current_label = 1;
coordList = zeros(8*size(K,1)*size(K,2),2);
cordnow = 1;
cordend = 1;


for i = 1:size(K,1)
    for j = 1:size(K,2)
        
        if(labels(i,j) == 0)
            coordList(cordend:(cordend+7),:) = [i-1 j-1; i-1 j; i-1 j+1;...
                         i j-1; i j+1; i+1 j-1; i+1 j; i+1 j+1];
            cordend = cordend+8;
            labels(i,j) = current_label;
            
            while(cordnow < cordend) 
                a = coordList(cordnow,1);
                b = coordList(cordnow,2);
                
                if(a>=1 && a<=size(K,1) && b>=1 && b<=size(K,2) ...
                        && labels(a,b) == 0 && K(a,b) == K(i,j))
                    labels(a,b) = current_label;
                    coordList(cordend:cordend+7,:) = [a-1 b-1; a-1 b; a-1 b+1;...
                                 a b-1; a b+1; a+1 b-1; a+1 b; a+1 b+1];
                    cordend = cordend+8;
                end
                
                cordnow = cordnow+1;
            end
            
            current_label = current_label + 1;
        end
    end
end


end

%% Old version that works fine but I'm trying to make it faster

% function [ labels] = connComp( K )
% %UNTITLED Summary of this function goes here
% %   Detailed explanation goes here
% 
% % Note: I use 8-way adjacency.
% 
% labels = zeros(size(K));
% current_label = 1;
% 
% for i = 1:size(K,1)
%     for j = 1:size(K,2)
%         
%         if(labels(i,j) == 0)
%             coordList = [i-1 j-1; i-1 j; i-1 j+1;...
%                          i j-1; i j+1; i+1 j-1; i+1 j; i+1 j+1];
%             labels(i,j) = current_label;
%             
%             while(size(coordList,1) > 1) 
%                 a = coordList(1,1);
%                 b = coordList(1,2);
%                 
%                 if(a>=1 && a<=size(K,1) && b>=1 && b<=size(K,2) ...
%                         && labels(a,b) == 0 && K(a,b) == K(i,j))
%                     labels(a,b) = current_label;
%                     coordList = [coordList; a-1 b-1; a-1 b; a-1 b+1;...
%                                  a b-1; a b+1; a+1 b-1; a+1 b; a+1 b+1];
%                 end
%                 
%                 coordList(1,:) = [];
%             end
%             
%             current_label = current_label + 1;
%         end
%     end
% end
% 
% 
% end
% 
%% Old version that I couldn't make work because of some bogus out of memory issue
% that really makes no sense because the internet is claiming that I'm
% passing by reference correctly (how do I check if I'm not?)

% current_label = 1;
% for i = 1:size(K,1)
%     for j = 1:size(K,2)
%         if(labels(i,j) == 0)
%             dfs(i,j,K,current_label,K(i,j),adjacency);
%             current_label = current_label + 1;
%         end
%     end
% end
% 
% toReturn = labels;
% 
% end
% 
% function [] = dfs(i,j,K,current_label,group,adjacency)
% global labels;
% 
% if(i >= 1 && j >= 1 && i <= size(K,1) && j <= size(K,2) ...
%         && labels(i,j) == 0 && group == K(i,j))
%     labels(i,j) = current_label;
%     
%     dfs(i+1,j,K,current_label,group,adjacency);
%     dfs(i,j+1,K,current_label,group,adjacency);
%     dfs(i-1,j,K,current_label,group,adjacency);
%     dfs(i,j-1,K,current_label,group,adjacency);
%     
%     if(adjacency == 8)
%         dfs(i-1,j-1,K,current_label,group,adjacency);
%         dfs(i-1,j+1,K,current_label,group,adjacency);
%         dfs(i+1,j-1,K,current_label,group,adjacency);
%         dfs(i+1,j+1,K,current_label,group,adjacency);
%     end
%     
% end
% 
% end