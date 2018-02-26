function [ newim ] = resizeIm( im, ratio, method )
% This is a really low level implementation of resize image
% I'm not doing anything clever here.
% In particular, it only works when ratio is a whole number bigger than 1.

% If the image sidelengths aren't perfectly divisible by the ratio, it
% chops off the far-side content

newim = zeros(floor(size(im,1)/ratio),floor(size(im,2)/ratio),size(im,3));

for i = 1:size(newim,1)
    for j = 1:size(newim,2)
        xlow = 1 + (i-1)*ratio;
        xhigh = i*ratio;
        ylow = 1 + (j-1)*ratio;
        yhigh = j*ratio;
        if(strcmp(method,'mean'))
            newim(i,j) = mean(reshape(im(xlow:xhigh,ylow:yhigh,:), [ratio*ratio size(im,3)]));
        else
            newim(i,j) = mode(reshape(im(xlow:xhigh,ylow:yhigh,:), [ratio*ratio size(im,3)]));
        end
    end
end

end

