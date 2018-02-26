% 
% imSize = [3,3];
% 
% im = zeros([imSize 4]);
% % im = randn(20,20,3)/10000;
% 
% for i = 1:imSize(1)
%     for j = 1:imSize(2)
%         if(i > imSize(1)/2)
%             im(i,j,1:3) = [1 0 0];
%         else
%             im(i,j,1:3) = [0 1 0];
%         end
%         if(j > imSize(2)/2)
%             im(i,j,4) = 1;
%         else
%             im(i,j,4) = 0;
%         end
%     end
% end
% 
% im = im + rand([imSize 4])/5;
% im(:,:,1:3) = rescaleIm(im(:,:,1:3));
% im(:,:,4) = rescaleIm(im(:,:,4));
% 
% subplot(1,2,1)
% imshow(im(:,:,1:3));
% subplot(1,2,2)
% imshow(im(:,:,4));
% imwrite(im,'~/Schoolwork/Chanussot/MBO_Code/data/test.tiff','tiff');

addpath ImageTools

A = load('~/Schoolwork/Chanussot/Matlab/SpectralClustering/Jadeplant/Data.mat');

im = rescaleIm(A.im);
im = nonlocal(im,4,3);
lidar = rescaleIm(A.lidar);
imwrite(cat(3,im,lidar),'~/Schoolwork/Chanussot/MBO_Code/data/jadeplant_both_nonlocal.tiff','tiff');

