function [signs, evecassign] = fixSigns(superU1,superU2)
% Input: superU1, superU2. Both are (# supervised points) x (# evecs)
%           matrices. Assumed to have rows registered correctly.
%     This function compares columns (aka eigenfunctions). Finds the
%     "optimal" arrangement (permutation and +- signs) to align
%     eigenvectors between sets.
% Output: evecassign: (# evecs) x 2 matrix. Each row is a pair of
%             eigenvectors (column 1 is set 1, column2 is set 2)
%         signs: corresponding +- adjustments for evecs

temp = superU2'*superU1;
tempsigns = sign(temp);
temp = temp .* tempsigns;
evecassign = Hungarian(temp);

signs = tempsigns(sub2ind(size(temp),evecassign(:,1),evecassign(:,2)));

% signs = zeros(1,size(U1,2));
% 
% for i = 1:size(super,1)
%     signs = signs + sign(U1(super(i,1),:) ./ U2(super(i,2),:));
% end
% 
% signs = sign(signs);

end
