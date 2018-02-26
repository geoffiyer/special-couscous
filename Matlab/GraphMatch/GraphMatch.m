%% I think this code is super old...

function [assign,Y1,Y2,Yassign] = GraphMatch(U1,U2,clusterSize)

%% Graph Laplacian and eigenvectors
% Eigenstuff using Nystrom 

%% Coarsify and match it up
% Still don't have a good choice of how to coarsify.
[Y1,Y1clusters] = CoarsifyGraph(U1,clusterSize,'rand');
[Y2,Y2clusters] = CoarsifyGraph(U2,clusterSize,'rand');

% Better way to do it: Select out only the parts of U we really need
smallU = U1(Y1clusters(:,1),:)*(U2(Y2clusters(:,1),:)');
Yassign = Hungarian(smallU);
assign = CoarseAssign(U1,U2,Y1clusters,Y2clusters,Yassign);

end
