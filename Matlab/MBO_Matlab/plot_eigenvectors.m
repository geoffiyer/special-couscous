
%% Plot some eigenvectors
% Load the eigenvectors from nystrom_result into evecs
% Plot them one by one (with a pause)
% Note: need to have image loaded in variable im becuase I use the size
%       to reshape

addpath '../ImageTools/'

if(~exist('im','var'))
    im = double(imread('~/Schoolwork/Chanussot/MBO_Code/data/umbrella_both.tiff'));
    lidar = rescaleIm(im(:,:,4));
    im = im(:,:,1:3);
end

evecs = dlmread('~/Schoolwork/Chanussot/MBO_Code/Output/nystrom_V.txt');

saveImages = true;
directory = '/home/gsiyer/Schoolwork/Chanussot/Papers/Image_Segmentation/Images/Umbrella/';

for i = 1:15
    tempmin = min(evecs(:,i));
    tempmax = max(evecs(:,i));
    temp = (evecs(:,i) - tempmin) / (tempmax - tempmin);
    if(i < 10)
        tempstr = strcat('evec0',num2str(i),'.png');
    else
        tempstr = strcat('evec',num2str(i),'.png');
    end
    imshow(64*reshape(temp, [size(im,2) size(im,1)])',jet);
    title(tempstr);
    
    if(saveImages)
        filename = strcat(directory,tempstr);
        imwrite(64*reshape(temp, [size(im,2) size(im,1)])', jet, filename, 'png');
    end
    
    pause;
end

close;
clear tempmin tempmax temp tempstr saveImages directory;