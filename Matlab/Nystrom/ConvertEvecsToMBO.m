
assert(logical(exist('I','var')),'Need eigenfunctions stored in variable I');
assert(logical(exist('D','var')),'Need eigenvalues stored in variable D');
assert(size(D,2) == 1,'Size of D is wrong!');
assert(size(D,1) >= 100, 'Need >= 100 eigenvalues to do the MBO');
assert(size(I,3) >= 100, 'Need >= 100 eigenfunctions to do the MBO');

tempI = reshape(permute(I,[2 1 3]), [size(I,1)*size(I,2) size(I,3)]);
tempI = tempI(:,1:100);

dlmwrite('~/Schoolwork/Chanussot/MBO_Code/Input/nystrom_V.txt',tempI,' ');
dlmwrite('~/Schoolwork/Chanussot/MBO_Code/Input/nystrom_D.txt',D,' ');