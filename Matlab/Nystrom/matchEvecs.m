function [V2] = matchEvecs(V1,V2)
% Input: V1, V2. Both are (# points) x (# evecs)
%           matrices. Assumed to have rows registered correctly.
%     This function compares columns (aka eigenfunctions). Finds the
%     "optimal" arrangement (permutation and +- signs) to align
%     eigenvectors between sets.
% Output: Edited version of V2 to have the proper match

temp = V2'*V1; % signed correllation matrix
tempsigns = sign(temp);
temp = temp .* tempsigns; % get rid of signs
evecassign = Hungarian(temp); % do the match

signs = tempsigns(sub2ind(size(temp),evecassign(:,1),evecassign(:,2))); % get relevant signs
V2 = V2.*signs'; % fix the +-
V2(:, evecassign(:,2)) = V2; % fix the ordering

end
