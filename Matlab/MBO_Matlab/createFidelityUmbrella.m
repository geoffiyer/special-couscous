
%% Fidelity for umbrella picture
% Note: this is about 5% fidelity

im = imread('~/Schoolwork/Chanussot/MBO_Code/data/umbrella_both.tiff');
imSize = size(im);

filename = '~/Schoolwork/Chanussot/Papers/Image_Segmentation/Images/Umbrella/fidelity.png';

fidelity = zeros([size(im,1) size(im,2)]);

         %x1  y1  x2  y2
boxes = [ 10  10  70  30; ...   white background wall
         200 150 230 180; ...   front umbrella
          50  70  80 100; ...   back umbrella
         320  20 350  50; ...   brown cabinet thingy
         100 200 120 220; ...   placeholder for dark-colored background objects
                          ...   (this is the blue chair at the bottom)
         170  10 190  30; ...   That gray bit of background wall near the top
        ];
    
for i = 1:size(boxes,1)
    fidelity(boxes(i,2):boxes(i,4), boxes(i,1):boxes(i,3)) = i/255;
end

imwrite(rescaleIm(fidelity)*64,jet,filename);

clear filename i boxes;
