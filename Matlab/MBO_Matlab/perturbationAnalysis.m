
%% Overall goal of this file is to do perturbation analysis using both spectral
% clustering and the MBO code. This is annoying because I have half C++ and
% half matlab.

%% NOTE: To make this work, I had to mess with libstdc++.so.6
% To fix, go to /usr/local/MATLAb/R2016b/sys/os/glnxa64/
% and rename libstdc++.so.6.old back to libstdc++.so.6

%% Make two different weight matrices, corresponding evecs

addpath ../SpectralClustering/
addpath ../ImageTools/
addpath ../Nystrom/

% A = load('DFC/DFC2013.mat');
% baseIm = A.im(:,:,[40,70,100]);
% nclasses = 6;

% A = load('Jadeplant/Data.mat');
% unionClasses1 = 3;
% unionClasses2 = 2;
% nclasses = 7;

A = load('../SpectralClustering/Umbrella/Data.mat');
nclasses = 6;  % nclasses hXas to match fidelity be careful
               % actually wait I think this numclasses is entirely ignored
               % The real numclasses is currently hardcoded into c++
               % Because I haven't bothered to mess with that yet
epsilon = 0.06;
numTrials = 1;

% A = load('DFC/DFC2015.mat');
% nclasses = 6;

p = inf;
im = A.im;
lidar = A.lidar;
[E,dex] = getTwoModalWeights(im, lidar, p);
[V,D] = getEvecs(E,dex,100);
% V = V./colnorm(V);

%% Poop first V to a text file, 
assert(size(D,2) == 1,'Size of D is wrong!');
assert(size(D,1) >= 100, 'Need >= 100 eigenvalues to do the MBO');
assert(size(V,2) >= 100, 'Need >= 100 eigenfunctions to do the MBO');

I = reshape(V,[size(im,1) size(im,2) size(V,2)]);
tempI = reshape(permute(I,[2 1 3]), [size(I,1)*size(I,2) size(I,3)]);
tempI = tempI(:,1:100);

dlmwrite('~/Schoolwork/Chanussot/MBO_Code/Input/From_Matlab/nystrom_V.txt',tempI,' ');
dlmwrite('~/Schoolwork/Chanussot/MBO_Code/Input/From_Matlab/nystrom_D.txt',D,' ');
imwrite(fidelity,'~/Schoolwork/Chanussot/MBO_Code/data/fidelityFromMatlab.tiff');

%% Run the c++ code on the first V. Save results in K
command = '/home/gsiyer/Schoolwork/Chanussot/MBO_Code/Geoff_Semisupervised_MBO/bin/a.exe';
[status,cmdout] = system(command);

imSize = [size(im,1) size(im,2) 4];
im = double(permute(reshape(dlmread('~/Schoolwork/Chanussot/MBO_Code/Output/image.txt'),[imSize(2),imSize(1),imSize(3)]),[2 1 3]));
lidar = rescaleIm(im(:,:,4));
im = rescaleIm(im(:,:,1:3));

K = reshape(dlmread('~/Schoolwork/Chanussot/MBO_Code/Output/classification_result.txt'), [size(im,2) size(im,1)])';

%% Now run the perturbation in a loop and compare.
KList = zeros(numTrials,1);

for i = 1:numTrials
    Epert = E + randn(size(E))*std(E(:))*epsilon;
    [Vpert,Dpert] = getEvecs(Epert,dex,100);
    % Vpert = Vpert./colnorm(Vpert);
    Vpert = matchEvecs(V,Vpert);
    
    assert(size(Dpert,2) == 1,'Size of D is wrong!');
    assert(size(Dpert,1) >= 100, 'Need >= 100 eigenvalues to do the MBO');
    assert(size(Vpert,2) >= 100, 'Need >= 100 eigenfunctions to do the MBO');
    
    Ipert = reshape(Vpert,[size(im,1) size(im,2) size(V,2)]);
    tempI = reshape(permute(Ipert,[2 1 3]), [size(I,1)*size(I,2) size(I,3)]);
    tempI = tempI(:,1:100);
    
    dlmwrite('~/Schoolwork/Chanussot/MBO_Code/Input/From_Matlab/nystrom_V.txt',tempI,' ');
    dlmwrite('~/Schoolwork/Chanussot/MBO_Code/Input/From_Matlab/nystrom_D.txt',Dpert,' ');
    
    %% Run the c++ code on Vpert. Save results in Kpert
    command = '/home/gsiyer/Schoolwork/Chanussot/MBO_Code/Geoff_Semisupervised_MBO/bin/a.exe';
    [status,cmdout] = system(command);
    
    imSize = [size(im,1) size(im,2) 4];
    im = double(permute(reshape(dlmread('~/Schoolwork/Chanussot/MBO_Code/Output/image.txt'),[imSize(2),imSize(1),imSize(3)]),[2 1 3]));
    lidar = rescaleIm(im(:,:,4));
    im = rescaleIm(im(:,:,1:3));
    fidelity = reshape(dlmread('~/Schoolwork/Chanussot/MBO_Code/Output/fidelity.txt'),[size(im,2) size(im,1)])';
    fidelity = rescaleIm(fidelity);
    
    Kpert = reshape(dlmread('~/Schoolwork/Chanussot/MBO_Code/Output/classification_result.txt'), [size(im,2) size(im,1)])';
    
    KList(i) = sum( K(:) ~= Kpert(:) ) / (size(im,1)*size(im,2));
end

%% plot some crap

subplot(1,3,1)
imshow(rescaleIm(K)*64,jet)
subplot(1,3,2)
imshow(rescaleIm(Kpert)*64,jet)
subplot(1,3,3)
imshow(K ~= Kpert)

temp = V./colnorm(V) - Vpert./colnorm(Vpert);
Vnorms = colnorm(temp);
figure
imshow(rescaleIm(reshape(V(:,1),[size(im,1) size(im,2)]))*64,jet);
figure
imshow(rescaleIm(reshape(Vpert(:,1),[size(im,1) size(im,2)]))*64,jet);

plotCrap = false;
if(plotCrap)
    if(exist('fidelity','var'))
        fidelity = rescaleIm(fidelity);
        lidar = rescaleIm(lidar);
        subplot(2,2,1);
        imshow(im);
        title('RGB image');
        subplot(2,2,2);
        imshow(lidar*64,jet);
        title('lidar image');
        subplot(2,2,3);
        imshow(fidelity*64,jet);
        title('fidelity');
        subplot(2,2,4);
        imshow(rescaleIm(K)*64,jet);
        title('classification');
    elseif(exist('lidar','var'))
        subplot(1,3,1);
        imshow(im);
        title('RGB image');
        subplot(1,3,2);
        imshow(lidar);
        title('lidar image');
        subplot(1,3,3);
        imagesc(K);
        title('classification');
    else
        subplot(1,2,1);
        imshow(im);
        title('RGB image');
        subplot(1,2,2);
        imagesc(K);
        title('Classification');
    end
end
