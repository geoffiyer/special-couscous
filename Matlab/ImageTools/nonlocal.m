
%% Nonlocal means
% Note on matlab reshape. It will do one column at a time.
% Example, A is 100 x 50 x 3
% I reshape it into B, 5000 x 3
% Then B(1:100,:) will be the same as A(:,1,:)
function [data] = nonlocal(data,squareRadius,sigma,imSize)

reshapeFlag = false;

if(size(data,3) > 1 || nargin < 4)
    imSize = size(data);
    if(size(imSize,2) == 2)
        imSize = [imSize 1];
    end
    reshapeFlag = true;
    data = reshape(data, [imSize(1)*imSize(2) imSize(3)]);
end

if(~exist('sigma','var'))
    sigma = 1;
end

tempData = data;

wd = 2*squareRadius + 1;

baseidx = zeros(wd*wd,1);
for i = 1:wd
    baseidx((1+(i-1)*wd):(i*wd)) = (1:wd) + (i-1)*imSize(1);
end

convvec = zeros(wd*wd, imSize(3));

for j = 1:wd
    for i = 1:wd
        convvec(i+(j-1)*wd,:) = repmat( exp( -norm( [i,j] - [(wd+1)/2 (wd+1)/2] ) / sigma) , [1 imSize(3)]);
    end
end

convvec = convvec / sum(abs(convvec(:,1)));

% (i,j) keeps track of the upper-left of the square
for j = (1-(wd-1)/2):(imSize(2) - (wd-1)/2)
    for i = (1-(wd-1)/2):(imSize(1) - (wd-1)/2)
        dataidx = baseidx + (i-1) + (j-1)*imSize(1);
        centeridx = dataidx((wd*wd+1)/2);
        % boundary cases
        if( i < 1 || j < 1 || i > (imSize(1)-wd+1) || j > (imSize(2)-wd+1) )
            newidx = false(wd,wd);
            for i2 = 0:(wd-1)
                for j2 = 0:(wd-1)
                    newidx(i2+1,j2+1) = (i + i2 > 0) && (j + j2 > 0) && (i + i2 <= imSize(1)) && (j + j2 <= imSize(2));
                end
            end
            newidx = reshape(newidx,[wd*wd 1]);
            data(centeridx,:) = sum( tempData(dataidx(newidx),:) .* convvec(newidx,:), 1);
            data(centeridx,:) = data(centeridx,:) / sum(convvec(newidx,1));
        else
            % pixel (i,j) in original image is represented as
            % index (j-1)*imSize(1) + i
            data(centeridx,:) = sum(tempData(dataidx,:) .* convvec, 1);
        end
    end
end

if(reshapeFlag)
    data = reshape(data, imSize);
end

end