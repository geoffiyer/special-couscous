
%% Fidelity for jadeplant picture
% roughly 5% fidelity

addpath ImageTools

im = double(imread('~/Schoolwork/Chanussot/MBO_Code/data/jadeplant_both.tiff'));
lidar = rescaleIm(im(:,:,4));
im = rescaleIm(im(:,:,1:3));
imSize = size(im);

filename = '~/Schoolwork/Chanussot/MBO_Code/data/fidelity_jadeplant_8class.tiff';

fidelity = zeros([size(im,1) size(im,2)]);

         %x1  y1  x2  y2
boxes = [ 30  30  50  50; ...   black background
         100 200 120 220; ...   front brown box
         100 260 120 280; ...   brown table
         300 170 320 190; ...   yellow cloth thing
         230 260 250 280; ...   block of wood under blue grid thingy
         200 200 220 220; ...   Pot for plant
                          ...     Note: this contains some of blue thing but that is fixed below
        ];
    
for k = 1:size(boxes,1)
    fidelity(boxes(k,2):boxes(k,4), boxes(k,1):boxes(k,3)) = k/255;
end

% Blue grid thing
idx = lidar > 0.9 & im(:,:,3) > im(:,:,2);
fidelity(idx) = (k+1)/255;

% Note: I decided to get rid of this and call the plant one whole object
% % Plant trunk
% idx = zeros(size(im,1),size(im,2));
% u = [181,237];        % start of line
% v = [134,260] - u;    % end of line, but centered at u now
% for i = 1:size(im,1)
%     for j = 1:size(im,2)
%         w = [i j] - u;
%         idx(i,j) = ((norm(w - dot(w,v)/dot(v,v)*v) < 3) && w(2) < v(2) && w(2) > 0);
%     end
% end
% idx = (idx == 1);
% fidelity(idx) = (k+2)/255;

% Plant
temp(1,1,1) = 0.6;
temp(1,1,2) = 0.6;
temp(1,1,3) = 0.1;
idx = (sqrt(sum( (im - temp).^2,3)) < 0.2 & lidar > 0.73 & lidar < 0.8);
idx(200:285,:) = 0;
fidelity(idx) = (k+2)/255;

% imshow(rescaleIm(fidelity)*64,jet);

imwrite(fidelity,filename);

clear filename i boxes idx idx2 temp;