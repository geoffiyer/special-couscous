
%% view MBO, and save pictures to file
% To run this you need the image loaded into the variable im
%  and obviously all the file input needs to work out

addpath ./ImageTools

if(~exist('im','var') || ~exist('lidar','var'))
    disp('Load im and lidar into matlab before running this function');
    return;
end

outDirectory = '~/Schoolwork/Chanussot/MeetingNotes/2017-08-10/MBO_Umbrella2Norm/';

% Load the eigenvectors from nystrom_result into evecs
% Save the first 100
% Note: need to have image loaded in variable im becuase I use the size
%       to reshape
evecs = dlmread('~/Schoolwork/Chanussot/MBO_Code/Output/nystrom_V.txt');

for i = 1:100
    tempmin = min(evecs(:,i));
    tempmax = max(evecs(:,i));
    temp = (evecs(:,i) - tempmin) / (tempmax - tempmin);
    if(i < 10)
        tempstr = strcat('evec0',num2str(i),'.png');
    else
        tempstr = strcat('evec',num2str(i),'.png');
    end
    
    filename = strcat(outDirectory,tempstr);
    imwrite(64*reshape(temp, [size(im,2) size(im,1)])', jet, filename, 'png');
    
end

clear tempmin tempmax temp tempstr i;

% rescale im and lidar
tempim = rescaleIm(im);
templidar = rescaleIm(lidar);

% Load fidelity into fidelity
fidelity = reshape(dlmread('~/Schoolwork/Chanussot/MBO_Code/Output/fidelity.txt'), [size(im,2) size(im,1)])';

% Load classification result into K
K = reshape(dlmread('~/Schoolwork/Chanussot/MBO_Code/Output/classification_result.txt'), [size(im,2) size(im,1)])';
K = smoothPixels(K,1);

% Sometimes images are scaled weird
% I prefer 0 to 1 (although jet uses 0 to 64)

% view classification
subplot(2,2,1);
imshow(tempim);
title('RGB image');
subplot(2,2,2);
imshow(templidar*64,jet);
title('lidar');
subplot(2,2,3);
imshow(rescaleIm(fidelity)*64,jet);
title('fidelity');
subplot(2,2,4);
imshow(K/max(max(K))*64,jet);
title('classification');

imwrite(tempim,strcat(outDirectory,'RGBimage.png'));
imwrite(rescaleIm(K)*64,jet,strcat(outDirectory,'classification.png'),'png');
imwrite(fidelity*64,jet,strcat(outDirectory,'fidelity.png'),'png');
if(exist('lidar','var'))
    imwrite(templidar*64,jet,strcat(outDirectory,'lidarimage.png'));
end
imwrite(borders(tempim,K),strcat(outDirectory,'imWithBorders.png'));
imwrite(cartoon(tempim,K),strcat(outDirectory,'cartoon.png'));

clear outDirectory filename tempim templidar;