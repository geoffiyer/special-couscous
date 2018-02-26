function [data] = zscore(data)

sig = std(data);
mu = mean(data);
data = (data - repmat(mu, [size(data,1) 1]))./repmat(sig, [size(data,1) 1]);

end