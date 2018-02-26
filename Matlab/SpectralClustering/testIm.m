%% Create a small test image so I can look at evecs

addpath ImageTools

% heightxwidth image. Try to make it as simple as possible
width = 10;
height = 10;
data = zeros(width,height,4);
for i=1:height
    for j=1:width
        if(i < 4)
            data(i,j,3) = 1;
        elseif (i > height-4)
            data(i,j,2) = 1;
        end
        if(i > height/2)
          data(i,j,4) = 1;
        end
    end
end
data = data + 0.01 * rand(width,height,4);
data = data - min(min(min(data))) + 0.01; % DON'T DIVIDE BY ZERO

% data1 = reshape(data(:,:,1:3),[9 3]);
% data2 = reshape(data(:,:,4),[9 1]);
% m = 9;
% p = inf;
% normType = 'angle';
% kernel.type = 'rbf';

data1 = data(:,:,1:3);
data2 = data(:,:,4);
nclasses = 4;
displayGraph = true;
pixelsmoothing = 0;
p = inf;

% Run the code on the dataset of choice
[error,I,K,bw,imWithBorders,scaling1,scaling2,E1,E2] = SpectralClassification(data1,data2,nclasses,displayGraph,pixelsmoothing,p);

for i=1:10
    subplot(2,5,i)
    imshow(rescaleIm(I(:,:,i))*64,jet)
end

clear displayGraph

%% No idea what this is

% Gray = A.lidar/255;
% RGB1 = cat(3, Gray, Gray, Gray);  % information stored in intensity
% RGB2 = Gray;
% RGB2(end, end, 3) = 0;  % All information in red channel
% GrayIndex = uint8(floor(Gray * 255));
% Map       = jet(255);
% RGB3      = ind2rgb(GrayIndex, Map);
% imshow(RGB3);
% 
% imwrite(RGB3,Map,'~/Schoolwork/Chanussot/Papers/ICIP_2017_Paper/Images/Umbrella/lidarColor.png');