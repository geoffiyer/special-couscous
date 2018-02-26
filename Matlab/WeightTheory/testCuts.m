
cuttype = 'ncut';
numPoints = 6;
data = [1  1  0 0;
        0 0.8 0.2 0;
        0 0 1 1;
        0.2 0 0 0.8];
     
W1 = weightMatrix(data,1);
WInf = weightMatrix(data,inf);

assign1 = bruteForceCut(W1,cuttype);
assignInf = bruteForceCut(WInf,cuttype);

if(size(data,2) == 2)
    subplot(1,2,1)
    hold on
    scatter(data(assign1,1),data(assign1,2),50,'blue','filled')
    scatter(data(~assign1,1),data(~assign1,2),50,'red','filled')
    title('1-norm segmentation')
    hold off
    subplot(1,2,2)
    hold on
    scatter(data(assignInf,1),data(assignInf,2),50,'blue','filled')
    scatter(data(~assignInf,1),data(~assignInf,2),50,'red','filled')
    title('Inf-norm segmentation')
    hold off
elseif(size(data,2) >= 3)
    subplot(1,2,1)
    hold on
    scatter3(data(assign1,1),data(assign1,2),data(assign1,3),50,'blue','filled')
    scatter3(data(~assign1,1),data(~assign1,2),data(~assign1,3),50,'red','filled')
    title('1-norm segmentation')
    hold off
    subplot(1,2,2)
    hold on
    scatter3(data(assignInf,1),data(assignInf,2),data(assignInf,3),50,'blue','filled')
    scatter3(data(~assignInf,1),data(~assignInf,2),data(~assignInf,3),50,'red','filled')
    title('Inf-norm segmentation')
    hold off
end
diffPoints = sum(abs(assign1 - assignInf));