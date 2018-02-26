
A = load('Umbrella/Data.mat');
nclasses = 8;
displayGraph = false;
pixelsmoothing = 2;
p = inf;

[~,I1,~,~,~,scaling1] = classify(A.im,A.im,nclasses,displayGraph,pixelsmoothing,p);
[~,I2,~,~,~,~,scaling2] = classify(A.lidar,A.lidar,nclasses,displayGraph,pixelsmoothing,p);

I = cat(3,I1(:,:,1:12),I2(:,:,1:12));
V = reshape(I, [size(I,1)*size(I,2) size(I,3)]);

maxIter = 100;
[idx, ~, ~ ] = eff_kmeans(V,nclasses,maxIter);
K = reshape(idx, [size(I,1), size(I,2)]);
K = smoothPixels(K,pixelsmoothing);
[bw,imWithBorders] = borders(A.im,K);
error = segmentationError(cat(3,A.im/scaling1,A.lidar/scaling2),K);  % Note: I messed around with dividing by dimension here

for i = [1 13 2 14 3 15 4  16 5 17]
    imshow(I(:,:,i)*64,jet);
    waitforbuttonpress;
end