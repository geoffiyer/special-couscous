%% After running spectral clustering, save results to a file
% To run this you should run the examples.m first and we can go from there

addpath ./ImageTools

im = A.im;
lidar = A.lidar;

outDirectory = '~/Schoolwork/Chanussot/MeetingNotes/2017-08-10/SpecClust_JadeplantNonlocal/';

% Leftover from the MBO version of this code. Doubt I'll ever use it but
% hey you never know
% ---------------------------------------
% % Load the eigenvectors from nystrom_result into evecs
% % Save the first 100
% % Note: need to have image loaded in variable im becuase I use the size
% %       to reshape
% evecs = dlmread('~/Schoolwork/Chanussot/MBO_Code/Output/nystrom_V.txt');
%
% for i = 1:100
%     tempmin = min(evecs(:,i));
%     tempmax = max(evecs(:,i));
%     temp = (evecs(:,i) - tempmin) / (tempmax - tempmin);
%     if(i < 10)
%         tempstr = strcat('evec0',num2str(i),'.png');
%     else
%         tempstr = strcat('evec',num2str(i),'.png');
%     end
%     
%     filename = strcat(outDirectory,tempstr);
%     imwrite(64*reshape(temp, [size(im,2) size(im,1)])', jet, filename, 'png');
%     
% end
% 
% clear tempmin tempmax temp tempstr i;

% ---------------------------------------

if(size(im,3)>3)
    lidar = im(:,:,4);
    im = im(:,:,1:3);
end

% Sometimes images are scaled weird
% I prefer 0 to 1 (although jet uses 0 to 64)
im = rescaleIm(im);
lidar = rescaleIm(lidar);

% view classification
if(exist('fidelity','var'))
    subplot(2,2,1);
    imshow(im);
    title('RGB image');
    subplot(2,2,2);
    imshow(lidar*64,jet);
    title('lidar');
    subplot(2,2,3)
    imshow(imWithBorders);
    subplot(2,2,4);
    imshow(K/max(max(K))*64,jet);
    title('classification');    
else
    subplot(1,3,1);
    imshow(im);
    title('RGB image');
    subplot(1,3,2);
    imshow(lidar*64,jet);
    title('lidar');
    subplot(1,3,3);
    imagesc(K);
    title('classification');
end

imwrite(im,strcat(outDirectory,'RGBimage.png'));
imwrite(K/max(max(K))*64,jet,strcat(outDirectory,'classification.png'),'png');
if(exist('lidar','var'))
    imwrite(lidar*64,jet,strcat(outDirectory,'lidarimage.png'));
end
imwrite(borders(im,K),strcat(outDirectory,'imWithBorders.png'));
imwrite(cartoon(im,K),strcat(outDirectory,'cartoon.png'));

clear outDirectory filename;
