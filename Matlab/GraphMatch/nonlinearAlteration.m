function [ data ] = nonlinearAlteration( data, picSize )
% I don't really know
%   mess with each data point in some nonlinear way
%   then see what happens to the graph match

dim = size(data,2);

data = data ./ max(data,[],1);

if(dim == 3)
    data = [cos(data(:,1)*pi/2) sin(data(:,2))*pi/2 data(:,3).^2];
end

xstart = round(picSize(2)*9/20);
xend   = round(picSize(2)*11/20);
ystart = round(picSize(1)*9/20);
yend   = round(picSize(1)*11/20);

end

