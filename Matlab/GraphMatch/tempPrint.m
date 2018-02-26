
foldername = '/home/gsiyer/Schoolwork/Chanussot/MeetingNotes/2017-11-20/';

pic1Size = size(pic1);
pic2Size = size(pic2);

pic1edit = [1 1 1];
pic2edit = [0.8 0.5 0.8];

temp = diffNorms(:,2) + diffNorms(:,3);
imwrite(reshape(pic1edit.*X1, pic1Size),strcat(foldername,'pictureX.png'));
imwrite(reshape(pic2edit.*X2, pic2Size),strcat(foldername,'pictureY.png'));
imwrite(reshape(pic2edit.*X2(assign1to2(:,2),:), pic2Size),strcat(foldername,'permY.png'));
imwrite(reshape(pic1edit.*X1(assign2to1(:,1),:), pic1Size),strcat(foldername,'permX.png'));
imwrite(reshape(rescaleIm(diffNorms(:,2))*64,[pic1Size(1) pic1Size(2)]),jet, ...
    strcat(foldername,'ChangesYtoX.png'));
imwrite(reshape(rescaleIm(diffNorms(:,3))*64,[pic1Size(1) pic1Size(2)]),jet, ...
    strcat(foldername,'ChangesXtoY.png'));
imwrite(reshape(rescaleIm(temp)*64,[pic1Size(1) pic1Size(2)]),jet, ...
    strcat(foldername,'combinedResult.png'));

clear temp foldername pic1edit pic2edit;