
addpath ../ImageTools

tic
displayGraph = true;
pixelsmoothing = 1;
p = inf;
saveResults = false;
directory = '~/Schoolwork/Chanussot/Papers/Image_Segmentation/Images/DFC2015/';

% Pick one of the following choices of dataset.

% % Jadeplant with nonlocal means
% A = load('Jadeplant/Data.mat');
% actualIm = A.im;
% A.im = nonlocal(A.im,3,2);
% nclasses = 9;

% % Umbrella:
% A = load('Umbrella/Data.mat');
% nclasses = 6;

% DFC2015:
A = load('DFC/DFC2015.mat');
actualIm = A.im;
A.im = nonlocal(A.im,3,2);
nclasses = 6;

% % Storage:
% A = load('Storage/Data.mat');
% nclasses = 8;

% Run the code on the dataset of choice
[error,I,K,scaling1,scaling2,E1,E2] = SpectralClassification(A.im,A.lidar,nclasses,displayGraph,pixelsmoothing,p);
toc
A.im = actualIm;
imWithBorders = borders(A.im,K);

if(saveResults)
    imwrite(rescaleIm(K)*64,jet,strcat(directory,'specClust.png'));
    imwrite(imWithBorders,strcat(directory,'specClustBorders.png'));
end

clear displayGraph pixelSmoothing p actualIm saveResults directory
