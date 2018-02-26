%% Look at evecs of facebook data graph.
%
% addpath ../GraphMatch
% edgeList = dlmread('./airport/USairport500.txt');
% edgeList = edgeList(:,1:2);
% edgeList = edgeList - min(min(edgeList)) + 1;
% edgeList = unique(sort(edgeList,2),'rows');
% G = graph(edgeList(:,1),edgeList(:,2));
%
% h1 = plot(G,'-ko','NodeLabel',[]);
% axis off;
% for i=1:size(h1.XData,2)
%    text(h1.XData(i)+0.1,h1.YData(i),num2str(i),'fontsize',32);
% end
% h1.MarkerSize = 12;
% h1.LineWidth = 4;
% h1.EdgeAlpha = 1;

% A = full(adjacency(G));
% sqrtD = sum(A,2).^(-1/2);
% D = sum(A,2);
% % L = eye(size(A)) - sqrtD.*A.*(sqrtD');
% L = diag(D)-A;
% [V,evals] = eigs(L,20,'SA');
% scatter(X(:,1),X(:,2),40,V(:,3),'filled');
%
% addpath ../Nystrom
% idx = eff_kmeans(V(:,1:5),5,200);
% fig = figure;
% % line([X(edgeList(:,1),1),X(edgeList(:,2),1)],[X(edgeList(:,1),2) X(edgeList(:,2),2)])
% hold on;
% for i=1:size(edgeList,1)
%     line([X(edgeList(i,1),1),X(edgeList(i,2),1)],[X(edgeList(i,1),2) X(edgeList(i,2),2)],'Color','black','LineWidth',1);
% end
% scatter(X(:,1),X(:,2), 80, idx, 'filled');
% hold off;
% axis off;
% fig = tightfig(fig);

%% Look at evecs of an unweighted graph.

filename = '/home/gsiyer/Schoolwork/Chanussot/Presentations/GISPA_2017-10-19/Images/SyntheticGraph/';

addpath ../GraphMatch
% circle1 = [cos((1:24)/24*2*pi)', sin((1:24)/24*2*pi)'];
% circle1 = [circle1; 0.8*cos((1:20)/20*2*pi)', 0.8*sin((1:20)/20*2*pi)'];
% circle2 = circle1 + [7 0];
% line = (1.5:0.5:5.5)';
% line = [line, zeros(size(line,1),1)];

% circle1 = randCircle([0 0], 2, 40);
% circle2 = randCircle([7 0], 2, 40);
% theline = randLine([2 0],[5,0], 1, 15);
% X = [circle1; theline; circle2];
clear circle1 circle2 theline

numNeighbors = 5;
nbds = kNearestNeighbors(X,X,numNeighbors+1);

idx1 = reshape(repmat(1:size(X,1),[numNeighbors 1]),[1 size(X,1)*numNeighbors]);
idx2 = reshape(nbds(:,(1:numNeighbors)+1)',[1 size(X,1)*numNeighbors]);
edgeList = [idx1',idx2'];
edgeList = unique(sort(edgeList,2),'rows');

G = graph(edgeList(:,1),edgeList(:,2));

h1 = plot(G,'-ko','Xdata',X(:,1),'Ydata',X(:,2),'NodeLabel',[]);
axis off;
% for i=1:size(h1.XData,2)
%    text(h1.XData(i)+0.1,h1.YData(i),num2str(i),'fontsize',32);
% end
h1.MarkerSize = 8;
h1.LineWidth = 1;
h1.EdgeAlpha = 1;
saveas(h1,strcat(filename,'graph.png'));

A = full(adjacency(G));
sqrtD = sum(A,2).^(-1/2);
D = sum(A,2);
% L = eye(size(A)) - sqrtD.*A.*(sqrtD');
L = diag(D)-A;
[V,evals] = eigs(L,20,'SA');
fig = figure;
hold on;
for i=1:size(edgeList,1)
    line([X(edgeList(i,1),1),X(edgeList(i,2),1)],[X(edgeList(i,1),2) X(edgeList(i,2),2)],'Color','black','LineWidth',1);
end
for j=1:20
    scatter(X(:,1),X(:,2),80,V(:,j),'filled');
    fig = tightfig(fig);
    axis off;
    if(j < 10)
        tempstr = strcat('evec0',num2str(j),'.png');
    else
        tempstr = strcat('evec',num2str(j),'.png');
    end
    saveas(fig,strcat(filename,tempstr));
end

addpath ../Nystrom
idx = eff_kmeans(V(:,1:4),4,200);
fig = figure;
hold on;
for i=1:size(edgeList,1)
    line([X(edgeList(i,1),1),X(edgeList(i,2),1)],[X(edgeList(i,1),2) X(edgeList(i,2),2)],'Color','black','LineWidth',1);
end
scatter(X(:,1),X(:,2), 80, idx, 'filled');
hold off;
axis off;
fig = tightfig(fig);
saveas(fig,strcat(filename,'segmentation.png'));

%% Make a stick figure in a 2d scatterplot.
% % Then look at laplacian eigenfunctions
%
% stickFigure =  [randLine([-1 -0.5],[0 1],0.2,200); ...
%     randLine([ 1 -0.5],[0 1],0.2,200); ...
%     randLine([ 0 1],[0 4],0.2,200); ...
%     randLine([-1 3],[1 3],0.2,200)];
% temp = rand(200,1)*2*pi;
% circle = 0.6*[cos(temp), sin(temp)] + [0 4.6]+0.2*(rand(200,2)-0.5);
% stickFigure = [stickFigure; circle];
% clear temp circle
%
% addpath ../Nystrom
% [V,L,scaling] = INys_SpectrEmbed(stickFigure, 200, 'L2');
% for i=1:10
%     scatter(stickFigure(:,1),stickFigure(:,2), 40, V(:,i), 'filled');
%     pause
% end
%
% idx = eff_kmeans(V(:,2:10),6,200);
% fig = figure;
% scatter(stickFigure(:,1),stickFigure(:,2), 40, idx, 'filled');
% axis off;
% fig = tightfig(fig);

%% Make a circle in a 2d scatterplot
% % Then look at laplacian eigenfunctions
%
% pointDensity = 2;
% numCircles = 100;
% disc = [0 0];
% for i = 1:numCircles
%     numPoints = round(i*2*pi*pointDensity);
%     disc = [disc; i*cos((1:numPoints)/numPoints*2*pi)', i*sin((1:numPoints)/numPoints*2*pi)'];
% end
% clear numPoints numCircles baseCircle i low high pointDensity
%
% addpath ../Nystrom
% [V,L,scaling] = INys_SpectrEmbed(disc, 200, 'L2');
% rmpath ../Nystrom
%
% filename = '/home/gsiyer/Schoolwork/Chanussot/Presentations/GISPA_2017-10-19/Images/DiscExample/ManifoldLaplacian/';
% fig = figure;
% colormap default;
% for i=2:3
%     scatter(disc(:,1),disc(:,2),100,V(:,i),'filled')
%     axis off;
%     fig = tightfig(fig);
%     if(i < 10)
%         tempstr = strcat('evec0',num2str(i),'.png');
%     else
%         tempstr = strcat('evec',num2str(i),'.png');
%     end
%     saveas(fig,strcat(filename,tempstr));
% end
%
% clear i filename tempstr
%
% theta = 175;
% theta = theta*pi/180;
% rot = [cos(theta) sin(theta); -sin(theta) cos(theta)];
% disc = disc * rot;

%% Two moons in a 2d scatterplot
%
% numPointsPerMoon = 400;
%
% num = numPointsPerMoon;
% xvals = (rand(num,1) - 0.5)*2;
% yvals = 1.3*(1 - xvals.^2)-0.3;
% noise = rand(num,2)*0.2;
% moon1 = [xvals yvals] + noise;
%
% xvals = (rand(num,1) - 0.5)*2;
% yvals = 1.3*(xvals.^2 - 1)+0.3;
% xvals = xvals + 1;
% noise = rand(num,2)*0.2;
% moon2 = [xvals yvals] + noise;
%
% moons = [moon1; moon2];
% clear xvals yvals noise moon1 moon2 numPointsPerMoon num
%
% addpath ../Nystrom
% [V,L,scaling] = INys_SpectrEmbed(moons, 200, 'L2');
%
% % scatter(moons(:,1),moons(:,2),80,V(:,2),'filled')
% [idx,~,~] = eff_kmeans(V(:,2:10),2,200);
% scatter(moons(:,1),moons(:,2),80,idx,'filled')