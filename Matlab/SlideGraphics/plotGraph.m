
% %% Example number 1
% X1 = [1 1; 0 0; 2 0; 0 2; 2 2];
% % X2 = [2 0; 1 1; 2 2; 0 2; 0 0];
% % super = [1 1];
%
% G = graph([1; 1; 1; 1; 2; 2; 3; 3; 4;],[2; 3; 4; 5; 3; 4; 4; 5; 5;]);
%
% h1 = plot(G,'Xdata',X1(:,1),'Ydata',X1(:,2),'NodeLabel',[]);
% axis off;
% for i=1:size(h1.XData,2)
%    text(h1.XData(i)+0.1,h1.YData(i),num2str(i),'fontsize',32);
% end
% h1.MarkerSize = 15;
% h1.LineWidth = 3;
%
% H = graph([2 2 2 2 1 5],[1 3 4 5 3 4]);
%
% h2 = plot(H,'Xdata',X2(:,1),'Ydata',X2(:,2),'NodeLabel',[]);
% axis off;
% for i=1:size(h2.XData,2)
%    text(h2.XData(i)+0.1,h2.YData(i),num2str(i),'fontsize',32);
% end
% h2.MarkerSize = 15;
% h2.LineWidth = 3;


%% Example number 2

filename = '/home/gsiyer/Schoolwork/Chanussot/Presentations/GISPA_2017-10-19/Images/GraphMatch/';
fontsize = 48;

G = graph([1 1 2 ], [2 3 3]);
G.Edges.Weight = [10 1 5]';

x = [0 1 0.5];
y = [0 0 0.7];

LWidths = 8*G.Edges.Weight/max(G.Edges.Weight);
figure;
h1 = plot(G,'-ko','Xdata',x,'Ydata',y,'Nodelabel',[],'LineWidth',LWidths);
for i=1:size(h1.XData,2)
    if(i==1)
        offset = -0.25;
    else
        offset = 0;
    end
    text(h1.XData(i)+0.1+offset,h1.YData(i),num2str(i),'fontsize',fontsize);
end
% text((h1.XData(1)+h1.XData(2))/2+0.0,...
%     (h1.YData(1)+h1.YData(2))/2+0.05,...
%     '10','fontsize',16);
h1.MarkerSize = 10;
h1.EdgeAlpha = 1;
axis off;
saveas(h1,strcat(filename,'isom1.png'));

H = graph([1 1 2 ], [2 3 3]);
H.Edges.Weight = [5 10 1]';

x = [0 0.5 0.1];
y = [0 0.1 0.3];

LWidths = 8*H.Edges.Weight/max(H.Edges.Weight);
figure;
h2 = plot(H,'-ko','Xdata',x,'Ydata',y,'Nodelabel',[],'LineWidth',LWidths);
textoffset = [-0.11 0; 0.05 0; 0.05 0.01];
for i=1:size(h2.XData,2)
    text(h2.XData(i)+textoffset(i,1),h2.YData(i)+textoffset(i,2),num2str(i),'fontsize',fontsize);
end
h2.MarkerSize = 10;
h2.EdgeAlpha = 1;
axis off;
saveas(h2,strcat(filename,'isom2.png'));
