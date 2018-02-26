%% Create data, set parameters

set1 = false;
set2 = false;
set3 = true;
GenerateData;
addpath ./general_routine
% parameters
N = 20;
r = 1;
NF = 100;
mu = 0.5;               %(1-mu)*L  + mu*(Ls)
options.graph.nn = 10;  %KNN graph number of neighbors
dim = 1;

r1 = []; rT1 = []; r2 = []; rT2 = [];
rl1 = []; rlT1 = []; rl2 = []; rlT2 = [];
rW1 = []; rWT1 = []; rW2 = []; rWT2 = [];
results = [];

%% Set up data (training and testing sets)

%50%-50% split for training and testing
XT1 = X1(1:2:end,:)';
YT1 = Y1(1:2:end,:);
T = length(XT1)/2;

Xtemp1 = X1(2:2:end,:);
Ytemp1 = Y1(2:2:end,:);

XT2 = X2(1:2:end,:)';
YT2 = Y2(1:2:end,:);

Xtemp2 = X2(2:2:end,:);
Ytemp2 = Y2(2:2:end,:);

[X1 Y1 U1 Y1U indices] = ppc(Xtemp1,Ytemp1,N,r);
[X2 Y2 U2 Y2U indices] = ppc(Xtemp2,Ytemp2,N,r);

X1 = X1';
X2 = X2';
U1 = U1(1:2:end,:)';
U2 = U2(1:2:end,:)';

clear *temp*

Y1U = zeros(length(U1),1);
Y2U = zeros(length(U2),1);

ncl = numel(unique(Y1));

Y = [Y1;Y1U;Y2;Y2U];
YT = [YT1;YT2];


[d1 n1] = size(X1);
[d2 n2] = size(X2);

[temp,u1] = size(U1);
[temp,u2] = size(U2);

n = n1+n2+u1+u2;
d = d1+d2;

n1=n1+u1;
n2=n2+u2;

[dT1 T1] = size(XT1);
[dT2 T2] = size(XT2);

dT = dT1+dT2;

%% Graph Laplacian Stuff

% 1) Data in a block diagonal matrix
Z = blkdiag([X1,U1],[X2,U2]); % (d1+d2) x (n1+n2)

% 2) graph Laplacians
G1 = buildKNNGraph([X1,U1]',options.graph.nn,1);
G2 = buildKNNGraph([X2,U2]',options.graph.nn,1);
W = blkdiag(G1,G2);
W = double(full(W));
clear G*

% Class Graph Laplacian
Ws = repmat(Y,1,length(Y)) == repmat(Y,1,length(Y))'; Ws(Y == 0,:) = 0; Ws(:,Y == 0) = 0; Ws = double(Ws);
Wd = repmat(Y,1,length(Y)) ~= repmat(Y,1,length(Y))'; Wd(Y == 0,:) = 0; Wd(:,Y == 0) = 0; Wd = double(Wd);


Sws = sum(sum(Ws));
Sw = sum(sum(W));
Ws = Ws/Sws*Sw;

Swd = sum(sum(Wd));
Wd = Wd/Swd*Sw;

Ds = sum(Ws,2); Ls = diag(Ds) - Ws;
Dd = sum(Wd,2); Ld = diag(Dd) - Wd;
D = sum(W,2); L = diag(D) - W;


% Tune the generalized eigenproblem
A = ((1-mu)*L  + mu*(Ls)); % (n1+n2) x (n1+n2) %
B = Ld;         % (n1+n2) x (n1+n2) %

%% KEMA - RBF KERNEL
disp('  Mapping with the RBF kernel ...')

% 2) Compute RBF kernels
sigma1 =  15*mean(pdist([X1]'));
K1 = kernelmatrix('rbf',[X1,U1],[X1,U1],sigma1);
sigma2 =  15*mean(pdist([X2]'));
K2 = kernelmatrix('rbf',[X2,U2],[X2,U2],sigma2);

K = blkdiag(K1,K2);

KT1 = kernelmatrix('rbf',[X1,U1],XT1,sigma1);
KT2 = kernelmatrix('rbf',[X2,U2],XT2,sigma2);


KAK = K*A*K;
KBK = K*B*K;

% 3) Extract all features (now we can extract n dimensions!)
[ALPHA LAMBDA] = gen_eig(KAK,KBK,'LM');

[LAMBDA j] = sort(diag(LAMBDA));
ALPHA = ALPHA(:,j);



% 3b) check which projections must be inverted (with the 'mean of projected
% samples per class' trick) and flip the axis that must be flipped
E1     = ALPHA(1:n1,:);
E2     = ALPHA(n1+1:end,:);
sourceXpInv = (E1'*K1*-1)';
sourceXp = (E1'*K1)';
targetXp = (E2'*K2)';


sourceXpInv = zscore(sourceXpInv);
sourceXp = zscore(sourceXp);
targetXp = zscore(targetXp);


ErrRec = zeros(numel(unique(Y1)),size(ALPHA,2));
ErrRecInv = zeros(numel(unique(Y1)),size(ALPHA,2));

m1 = zeros(numel(unique(Y1)),size(ALPHA,2));
m1inv = zeros(numel(unique(Y1)),size(ALPHA,2));
m2 = zeros(numel(unique(Y1)),size(ALPHA,2));

cls = unique(Y1);

for j = 1:size(ALPHA,2)
    
    for i = 1:numel(unique(Y1))
        
        m1inv(i,j) = mean(sourceXpInv([Y1;Y1U]==cls(i),j));
        m1(i,j) = mean(sourceXp([Y1;Y1U]==cls(i),j));
        m2(i,j) = mean(targetXp([Y2;Y2U]==cls(i),j));
        
        ErrRec(i,j) = sqrt((mean(sourceXp([Y1;Y1U]==cls(i),j))-mean(targetXp([Y2;Y2U]==cls(i),j))).^2);
        ErrRecInv(i,j) = sqrt((mean(sourceXpInv([Y1;Y1U]==cls(i),j))-mean(targetXp([Y2;Y2U]==cls(i),j))).^2);
        
    end
end


mean(ErrRec);
mean(ErrRecInv);

Sc = max(ErrRec)>max(ErrRecInv);
ALPHA(1:n1,Sc) = ALPHA(1:n1,Sc)*-1;

% 4) Project the data
nVectRBF = min(NF,rank(KBK));
nVectRBF =  min(nVectRBF,rank(KAK));

for Nf = 1:nVectRBF
    
    E1     = ALPHA(1:n1,1:Nf);
    E2     = ALPHA(n1+1:end,1:Nf);
    Phi1toF = E1'*K1;
    Phi2toF = E2'*K2;
    
    Phi1TtoF = E1'*KT1;
    Phi2TtoF = E2'*KT2;
    
    % 5) IMPORTAT: Normalize!!!!
    m1 = mean(Phi1toF');
    m2 = mean(Phi2toF');
    s1 = std(Phi1toF');
    s2 = std(Phi2toF');
    
    Phi1toF = zscore(Phi1toF')';
    Phi2toF = zscore(Phi2toF')';
    
    Phi1TtoF = ((Phi1TtoF' - repmat(m1,2*T,1))./ repmat(s1,2*T,1))';
    Phi2TtoF = ((Phi2TtoF' - repmat(m2,2*T,1))./ repmat(s2,2*T,1))';
    
    
    
    
    % 6) Predict
    Ypred           = classify([Phi1toF(:,1:ncl*N)]',[Phi1toF(:,1:ncl*N),Phi2toF(:,1:ncl*N)]',[Y1;Y2]);
    Reslatent1Kernel2 = assessment(Y1,Ypred,'class');
    
    Ypred           = classify([Phi1TtoF]',[Phi1toF(:,1:ncl*N),Phi2toF(:,1:ncl*N)]',[Y1;Y2]);
    Reslatent1Kernel2T = assessment(YT1,Ypred,'class');
    
    Ypred           = classify([Phi2toF(:,1:ncl*N)]',[Phi1toF(:,1:ncl*N),Phi2toF(:,1:ncl*N)]',[Y1;Y2]);
    Reslatent2Kernel2 = assessment(Y2,Ypred,'class');
    
    Ypred           = classify([Phi2TtoF]',[Phi1toF(:,1:ncl*N),Phi2toF(:,1:ncl*N)]',[Y1;Y2]);
    Reslatent2Kernel2T = assessment(YT2,Ypred,'class');
    
    
    r1 = [r1; Reslatent1Kernel2.OA];
    rT1 = [rT1; Reslatent1Kernel2T.OA];
    
    r2 = [r2; Reslatent2Kernel2.OA];
    rT2 = [rT2; Reslatent2Kernel2T.OA];
    
end

results.RBF{r,dim}.X1 = r1;
results.RBF{r,dim}.XT1 = rT1;
results.RBF{r,dim}.X2 = r2;
results.RBF{r,dim}.XT2 = rT2;


%% unprojected

%train error
YpredO11 = classify(X1',X1',Y1);
ResOrig11 = assessment(Y1,YpredO11,'class');
r_la11 = ResOrig11.OA;

YpredO22 = classify(X2',X2',Y2);
ResOrig22 = assessment(Y2,YpredO22,'class');
r_la22 = ResOrig22.OA;

%test error
YpredT11 = classify(XT1',X1',Y1);
ResT11 = assessment(YT1,YpredT11,'class');
r_un11 = ResT11.OA;

YpredT22 = classify(XT2',X2',Y2);
ResT22 = assessment(YT2,YpredT22,'class');
r_un22 = ResT22.OA;


results.Upper{r,dim}.X1= r_la11;
results.Upper{r,dim}.XT1 = r_un11;
results.Upper{r,dim}.X2 = r_la22;
results.Upper{r,dim}.XT2 = r_un22;


if size(XT1,1) < size(XT2,1)
    XT1 = [XT1; 0.5+zeros(1,length(XT1))];
end

if size(XT2,1) < size(XT1,1)
    XT2 = [XT2; 0.5+zeros(1,length(XT2))];
end

% PLOT 1: original data
if min(size(XT1,1),size(XT2,1)) == 2
    
    figure,
    subplot(1,2,1)
    scatter(XT1(1,:),XT1(2,:),20,YT1,'f'), hold on, scatter(XT2(1,:),XT2(2,:),20,YT2),colormap(jet)
    title('original data (colors are classes)')
    grid on
    
    
    subplot(1,2,2)
    plot(XT1(1,:),XT1(2,:),'r.'), hold on, plot(XT2(1,:),XT2(2,:),'.'),colormap(jet)
    %legend('Domain 1','Domain 2')
    grid on
    title('Domains (red = X1, blue= X2)')
    
    
else
    
    figure,
    subplot(1,2,1)
    scatter3(XT1(1,:),XT1(2,:),XT1(3,:),20,YT1,'f'), hold on, scatter3(XT2(1,:),XT2(2,:),XT2(3,:),20,YT2),colormap(jet)
    title('original data (colors are classes)')
    grid on
    axis image
    
    subplot(1,2,2)
    plot3(XT1(1,:),XT1(2,:),XT1(3,:),'r.'), hold on, plot3(XT2(1,:),XT2(2,:),XT2(3,:),'.'),colormap(jet)
    %legend('Domain 1','Domain 2')
    grid on
    title('Domains (red = X1, blue= X2)')
    axis image
    
end

figure
subplot(1,2,1)
scatter(Phi1TtoF(1,:),Phi1TtoF(2,:),20,YT1,'f'), hold on, scatter(Phi2TtoF(1,:),Phi2TtoF(2,:),20,YT2),colormap(jet),hold off
title('Projected data (RBF)'),grid on
axis([-2.5 2.5 -2.5 2.5])

subplot(1,2,2)
scatter(Phi1TtoF(1,:),Phi1TtoF(2,:),45,'r+'), hold on, scatter(Phi2TtoF(1,:),Phi2TtoF(2,:),45,'bo'),colormap(jet),hold off
title('Projected data (RBF, domains)'),grid on
axis([-2.5 2.5 -2.5 2.5])