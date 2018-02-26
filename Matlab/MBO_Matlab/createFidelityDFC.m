%% Fidelity for DFC picture

if(~exist('im','var'))
    disp('Create an im before making fidelity');
else

filename = '~/Schoolwork/Chanussot/MBO_Code/data/fidelity_DFC_6classSMOOTH.tiff';

fidelity = zeros([size(im,1) size(im,2)]);

addpath ../ImageTools
lidar = rescaleIm(lidar);
im = rescaleIm(im);

% Tops of roofs
idx = (lidar > 0.93);
fidelity(idx) = 1/255;

% Mid-level roofs
idx = lidar > 0.295 & lidar < 0.31;
fidelity(idx) = 2/255;
% box: x1 y1 x2 y2
% box = [30 300 60 330];
% fidelity(box(2):box(4),box(1):box(3)) = 2/255;
% clear box;

% Grass
% box = [86 165 105 180];
% fidelity(box(2):box(4),box(1):box(3)) = 3/255;
idx = (im(:,:,2) > im(:,:,3)+0.12) & (abs(im(:,:,2) - im(:,:,1)) < 0.030);
fidelity(idx) = 3/255;
clear box;

% Light Road
idx = zeros(size(im,1),size(im,2));
u = [1,344];        % start of line
v = [204,394] - u;  % end of line, but centered at u now
for i = 1:size(im,1)
    for j = 1:size(im,2)
        w = [i j] - u;
        idx(i,j) = (norm(w - dot(w,v)/dot(v,v)*v) < 3);
    end
end
idx = logical(idx);
fidelity(idx) = 4/255;

% Dark Road
idx = zeros(size(im,1),size(im,2));
u = [250,90];        % start of line
v = [191,275] - u;  % end of line, but centered at u now
for i = 1:size(im,1)
    for j = 1:size(im,2)
        w = [i j] - u;
        idx(i,j) = ((norm(w - dot(w,v)/dot(v,v)*v) < 3) && w(2) < v(2) && w(2) > 0);
    end
end
idx = logical(idx);
fidelity(idx) = 5/255;

% White part
idx = (mean(im,3) > 0.8);
fidelity(idx) = 6/255;

fidpercent = size(fidelity(fidelity~=0))/(imSize(1)*imSize(2));
fidpercent = fidpercent(1);

fidelity = smoothPixels(fidelity,3);

imshow(rescaleIm(fidelity)*64,jet)
imwrite(fidelity,filename);

clear filename i j boxes idx idx2 temp imnorm u v w;

end