

% Small test image
im = double(imread('~/Schoolwork/Chanussot/Geoff_MBO_Code/data/test.tiff'));
lidar = rescaleIm(im(:,:,4));
im = im(:,:,1:3);

% umbrella
addpath('../ImageTools');
im = imread('~/Schoolwork/Chanussot/MBO_Code/data/umbrella_both.tiff');
imSize = size(im);
im = double(permute(reshape(dlmread('~/Schoolwork/Chanussot/MBO_Code/Output/image.txt'),[imSize(2),imSize(1),imSize(3)]),[2 1 3]));
lidar = im(:,:,4);
im = im(:,:,1:3);
fidelity = reshape(dlmread('~/Schoolwork/Chanussot/MBO_Code/Output/fidelity.txt'),[size(im,2) size(im,1)])';
fidelity = rescaleIm(fidelity);

% cow
im = double(imread('~/Schoolwork/Chanussot/Geoff_MBO_Code/data/cow.tiff'));

% DFC
im = double(imread('~/Schoolwork/Chanussot/MBO_Code/data/DFC_both.tiff'));
imSize = size(im);
im = double(permute(reshape(dlmread('~/Schoolwork/Chanussot/MBO_Code/Output/image.txt'),[imSize(2),imSize(1),imSize(3)]),[2 1 3]));
lidar = im(:,:,4);
im = im(:,:,1:3);
fidelity = reshape(dlmread('~/Schoolwork/Chanussot/MBO_Code/Output/fidelity.txt'),[size(im,2) size(im,1)])';

% jade plant
addpath('./ImageTools');
im = imread('~/Schoolwork/Chanussot/MBO_Code/data/jadeplant_both.tiff');
imSize = size(im);
im = double(permute(reshape(dlmread('~/Schoolwork/Chanussot/MBO_Code/Output/image.txt'),[imSize(2),imSize(1),imSize(3)]),[2 1 3]));
lidar = im(:,:,4);
im = im(:,:,1:3);
fidelity = reshape(dlmread('~/Schoolwork/Chanussot/MBO_Code/Output/fidelity.txt'),[size(im,2) size(im,1)])';
fidelity = rescaleIm(fidelity);