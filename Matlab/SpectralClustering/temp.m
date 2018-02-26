figure
for i=1:15
    subplot(3,5,i);
    imshow(I(:,:,i)*64,jet);
end