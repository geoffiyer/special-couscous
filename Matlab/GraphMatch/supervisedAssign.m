function [ assign ] = supervisedAssign( U, neighbors1, neighbors2, super )
% Match U2 to U1 using some stuff I'll write a better comment later

% Currently each node is counted as its own neighbor. I guess I'm okay with
% that and will keep it for now.

assign = [];

% Go through each supervised match and assign 5 neighbors.
% For now we are okay with many-to-many matchings
for i = 1:size(super,1)
    Utemp = U(neighbors2(super(i,2),:),neighbors1(super(i,1),:));
    inv_assign = munkres(1-Utemp);
    % temp_assign = [inv_assign', (1:size(inv_assign,2))'];  % not needed
    assign = [assign; neighbors1(super(i,1),inv_assign)', neighbors2(super(i,2),:)'];
end


end