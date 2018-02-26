%% NOTE: This doesn't do anything right now. Just a space filler for if I
% decide to really implement this dude.

function [ assign ] = HierarchicalAssign(U1, U2)
%HierarchicalAssign Assign at coarse level first, then at fine level.
%   Input: U1 = eigenvectors from X1
%          U2 = eigenvectors from X2
%             (compare these to create edgeweights)
%   Output: Assignment between X1 and X2



clustSize = 50;
% Still don't have a good choice of how to coarsify.
[~,Y1clusters] = CoarsifyGraph(U1,clustSize,'rand');
[~,Y2clusters] = CoarsifyGraph(U2,clustSize,'rand');

% Better way to do it: Select out only the parts of U we really need
Usmall = real(U1(Y1clusters(:,1),:) * U2(Y2clusters(:,1),:)');
Yassign = Hungarian(Usmall);
assign = CoarseAssign(U1,U2,Y1clusters,Y2clusters,Yassign);

end

