
A = dlmread('../../Graph_PATH_Code/graphm/test_Geoff/m_a_1EWK');
B = dlmread('../../Graph_PATH_Code/graphm/test_Geoff/m_a_1U19');
super = [1 14; 2 2];

degsA = sum(A,2);
A = eye(size(A)) - (1./sqrt(degsA)).*A.*(1./sqrt(degsA))';

degsB = sum(B,2);
B = eye(size(B)) - (1./sqrt(degsB)).*B.*(1./sqrt(degsB))';

[U1,D1] = eig(A,'vector');
[U2,D2] = eig(B,'vector');
U1 = U1(:,2:6);
U2 = U2(:,2:6);
[signs, evecassign] = fixSigns(U1(super(:,1),:),U2(super(:,2),:));
U2 = U2.*signs';
U2(:, evecassign(:,2)) = U2;

[assign, Y1, Y2, Yassign] = GraphMatch(U1,U2,1);

normSquareError = matchNormSquareError(assign,A,B);

% viewMatch(assign1to2,X1,X2);