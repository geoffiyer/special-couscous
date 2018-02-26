function X = load_data(filedir,filename, patch)
%% filedir is what it sounds like
% filename is without the extension for some silly reason
% patch is in the format [x1 x2 y1 y2]

%% Do all the work
imgpath = strcat(pwd,'/grss_dfc_2009/',filedir,filename,'.img');
hdrpath = strcat(pwd,'/grss_dfc_2009/',filedir,filename,'.hdr');

addpath('envi');
[X,~] = enviread(imgpath,hdrpath);

rmpath('envi');

if(nargin < 3)
    patch = [1 size(X,2) 1 size(X,1)];
end

xmin = max(patch(1),1);
xmax = min(patch(2), size(X,2));
ymin = max(patch(3),1);
ymax = min(patch(4), size(X,1));

X = X(ymin:ymax, xmin:xmax, :);

end