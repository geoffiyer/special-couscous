%% Outdated.
% I just hardcoded it for this specific dataset.
% Screw trying to be smart.

function [ im ] = killOutliers( im )
% Super basic outlier-removal function
% Nothing fancy, no defining what "outlier" really means.
% This is for just some testing really.
% Idea: Reduce the max by setting all points above lamda to lamda

addpath ../ImageTools/

if(nargin < 2)
    recursive = true;
end
im = rescaleIm(double(im));

while(recursive && (median(im(:)) < 0.02 || median(im(:)) > 0.98) )
    classes = otsu(im,2);
    
    oneidx = find(classes == 1,1);
    twoidx = find(classes == 2,1);
    if(im(oneidx) > im(twoidx))
        classes = 3-classes;
    end
    
    oneSize = size(classes(classes == 1),1);
    twoSize = size(classes,1)*size(classes,2) - oneSize;
    
    if oneSize > twoSize
        im(classes == 2) = min(im(classes == 2));
    else
        im(classes == 1) = max(im(classes == 1));
    end
    im = rescaleIm(double(im));
end

% im = rescaleIm(double(im));
% im = -2*im.^3 + 3*im.^2;
% 
% im = rescaleIm(double(im));
% im = 4*(im-1/2).^3 + 1/2;

end

