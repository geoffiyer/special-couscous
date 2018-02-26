
GenerateData;
X1 = X1(1:2,:);
% X_1 = X_1(:,[1 3 2 4 6 5 7 9 8]);
X2 = X2(1:2,:);

% % Trying an easier data set to understand better.
% X1 = [1 1.1 0.9 0 0 0; 0 0 0 1 1.1 0.9];
% X2 = [0 1 0 0.9 0 1.1; 1 0 1.1 0 0.9 0];

[nrx1, ncx1] = size(X1);
[nrx2, ncx2] = size(X2);
d = 2;   % number of eigenvectors to use
         % can change however I want
d = min([d, nrx1, nrx2]);

% I coded this myself, but it's buggy:
% S_ss = X_1*X_1';
% S_st = X_1*X_2';
% S_tt = X_2*X_2';
% 
% A = (S_ss)^-1*S_st*(S_tt)^-1*S_st';
% 
% [V,D] = eig(A);
% u_s = V;
% u_t = (S_tt)^-1*S_st*u_s*D^-1;
% u_t = u_t./repmat(sqrt(sum(u_t.^2,1)),numRows,1);

[u_s, u_t, r] = cca(X1,X2);
[r,I] = sort(r,'descend');

P1 = zeros(d,ncx1);
P2 = zeros(d,ncx2);

for i = 1:d
    P1 = P1 + dot(repmat(u_s(:,I(i)),1,ncx1),X1).*u_s(:,I(i));
    P2 = P2 + dot(repmat(u_t(:,I(i)),1,ncx2),X2).*u_t(:,I(i));
end

figure
hold on
myscat(P1,'o','blue');
myscat(P2,'+','red');
legend('P1','P2','Location','bestoutside');
hold off

ax = xlim;
ay = ylim;
s = max(ax(2)-ax(1),ay(2)-ay(1))/2;
xlim([(ax(1)+ax(2))/2 - s, (ax(1)+ax(2))/2 + s]);
ylim([(ay(1)+ay(2))/2 - s, (ay(1)+ay(2))/2 + s]);

% scale [a,b] to [A,B]
% (a+b)/2 remains fixed.


% 11/10, 12pm: My code is still buggy :(. I grabbed a working version of
% CCA from the internet and am using that now. Order of the points
% definitely matters.
%   The good news is that with the extra assumption of order, CCA along
% with residue idea is doing a good job of picking out the "outlier"
% points. Can look at X_1 - proju_s and see which dudes were not correctly
% classified, then maybe do something with this.

% 11/10, 1am: Found a pretty important bug and fixed it (scaled u_t wrong). Below comments are
% invalidated. Tried again. Still had similar behavior though. Both X_1 and
% X_2 are not represented well by the projection. I'm going to bed now
% For tomorrow: See to what extent the ordering of the points matters. It
% could be that this entire correlation thing is a "waste of time" if it
% assumed I've already made a bijection.