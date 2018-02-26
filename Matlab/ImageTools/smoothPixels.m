function [K2] = smoothPixels(K,squaresize)

n = squaresize;
K2 = K;

if(squaresize ~= 0)
    
    for i = 1:size(K,1)
        for j = 1:size(K,2)
            xmin = max(1,i-n);
            xmax = min(size(K,1),i+n);
            ymin = max(1,j-n);
            ymax = min(size(K,2),j+n);
            K2(i,j) = mode(reshape(K(xmin:xmax, ymin:ymax),[(xmax-xmin+1)*(ymax-ymin+1) 1]));
        end
    end
end

end