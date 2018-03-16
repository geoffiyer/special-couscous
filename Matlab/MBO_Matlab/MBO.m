function [ K, info ] = MBO(V, D, numClasses, fidelity, fidelityConst, dtConst)
%MBO Summary of this function goes here
%   Detailed explanation goes here

if(nargin < 5)
    fidelityConst = 1;
end

if(nargin < 6)
    dtConst = 1;
end

assert(size(D,2) == 1,'Size of D is wrong!');
assert(size(D,1) == size(V,2),'Size of D does not match size of V')

I = reshape(V,[size(fidelity,1) size(fidelity,2) size(V,2)]);
I = reshape(permute(I,[2 1 3]), [size(I,1)*size(I,2) size(I,3)]);
L = [numClasses; size(V,2); fidelityConst; dtConst];

dlmwrite('~/Schoolwork/Chanussot/MBO_Code/Input/From_Matlab/nm.txt',L,' ');
dlmwrite('~/Schoolwork/Chanussot/MBO_Code/Input/From_Matlab/nystrom_V.txt',I,' ');
dlmwrite('~/Schoolwork/Chanussot/MBO_Code/Input/From_Matlab/nystrom_D.txt',D,' ');
if(max(fidelity(:))*255 ~= numClasses)
    warning('fidelity does not match numClasses. Maybe you forgot to rescale??')
end
imwrite(fidelity,'~/Schoolwork/Chanussot/MBO_Code/data/fidelityFromMatlab.tiff');

command = '/home/gsiyer/Schoolwork/Chanussot/MBO_Code/Geoff_Semisupervised_MBO/bin/a.exe';
[status,info] = system(command);

try
    K = reshape(dlmread('~/Schoolwork/Chanussot/MBO_Code/Output/classification_result.txt'), [size(fidelity,2) size(fidelity,1)])';
catch
    K = zeros(size(fidelity));
end

end

