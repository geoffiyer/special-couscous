

x1 = 0;
xmargin = 5;
x1len = 100;
x2len = 50;
x3len = 150;
x4len = 50;
x5len = 150;
xtotal = x1 + x1len + x2len + x3len + x4len + x5len + 4*xmargin;

yint = 30;
ymargin = 10;

dotsize = 2;

arrowmargin = 8;

fontsize = 14;

y5 = 0;
y4 = y5 + yint + ymargin;
y3 = y4 + yint + ymargin;
y2 = y3 + yint + ymargin;
y1 = y2 + yint + ymargin;
ytotal = y1 + yint;

x2 = x1+x1len+xmargin;
x3 = x2+x2len+xmargin;
x4 = x3+x3len+xmargin;
x5 = x4+x4len+xmargin;

h = figure;
axis off
axis([x1 xtotal y5 ytotal]);

rectangle('Position',[x1 y1 x1len yint], 'Curvature', [0.1 0.3]);
text(x1+4.5*xmargin,y1+yint/2,'Data 1','FontSize',fontsize);
rectangle('Position',[x1 y2 x1len yint], 'Curvature', [0.1 0.3]);
text(x1+4.5*xmargin,y2+yint/2,'Data 2','FontSize',fontsize);
rectangle('Position',[x1 y3 x1len yint], 'Curvature', [0.1 0.3]);
text(x1+4.5*xmargin,y3+yint/2,'Data 3','FontSize',fontsize);
rectangle('Position',[x1 y5 x1len yint], 'Curvature', [0.1 0.3]);
text(x1+4.5*xmargin,y5+yint/2,'Data n','FontSize',fontsize);

dot1y = y4;
dot2y = y4 + yint/2 - dotsize/2;
dot3y = y3 - ymargin - dotsize;

rectangle('Position',[ x1+x1len/2 dot1y dotsize dotsize], 'FaceColor', 'black');
rectangle('Position',[ x1+x1len/2 dot2y dotsize dotsize], 'FaceColor', 'black');
rectangle('Position',[ x1+x1len/2 dot3y dotsize dotsize], 'FaceColor', 'black');

arrow([x2 y1+yint/2],[x2+x2len y3+yint/2 + 2*arrowmargin]);
arrow([x2 y2+yint/2],[x2+x2len y3+yint/2 +   arrowmargin]);
arrow([x2 y3+yint/2],[x2+x2len y3+yint/2]                );
arrow([x2 y5+yint/2],[x2+x2len y3+yint/2 -   arrowmargin]);
text(x2+x2len/2, y4, 'Data Fusion', 'Fontsize', fontsize);

rectangle('Position',[x3 (y4+yint/2) x3len (2*yint+2*ymargin)], 'Curvature', [0.2 0.3]);
text(x3+3*xmargin,y3+yint/2,'Feature Space','FontSize',fontsize);

arrow([x4 y3+yint/2],[x4+x4len y3+yint/2]);
arrow([x4+x4len/2 y4+ymargin/2], [x4+x4len/2 y3+yint/2 - ymargin/2]);
text(x4+x4len/2, y4, 'Further Processing','FontSize',fontsize);

rectangle('Position',[x5 y4+yint/2+ymargin x5len 2*yint], 'Curvature', [0.2 0.3]);
text(x5+5*xmargin,y3+yint/2,'Final Result','FontSize',fontsize);
