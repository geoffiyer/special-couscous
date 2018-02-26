function im = borders(im,K)

borderMat = zeros(size(im,1),size(im,2));

for i = 2:(size(K,1)-1)
    for j = 2:(size(K,2)-1)
        
        if( max( [abs(K(i,j) - K(i-1,j)), ...
                abs(K(i,j) - K(i+1,j)), ...
                abs(K(i,j) - K(i,j-1)), ...
                abs(K(i,j) - K(i,j+1))] )  > 0 )
            borderMat(i,j) = 1;
            im(i,j,:) = [1, zeros(1,size(im,3)-1)];
        end
    end
end

i = 1;
for j = 2:(size(K,2)-1)
    if( max( [ abs(K(i,j) - K(i+1,j)), ...
            abs(K(i,j) - K(i,j-1)), ...
            abs(K(i,j) - K(i,j+1))] )  > 0 )
        borderMat(i,j) = 1;
        im(i,j,:) = [1, zeros(1,size(im,3)-1)];
    end
end

i = size(im,1);
for j = 2:(size(K,2)-1)
    
    if( max( [abs(K(i,j) - K(i-1,j)), ...
            abs(K(i,j) - K(i,j-1)), ...
            abs(K(i,j) - K(i,j+1))] )  > 0 )
        borderMat(i,j) = 1;
        im(i,j,:) = [1, zeros(1,size(im,3)-1)];
    end
end

j = 1;
for i = 2:(size(K,1)-1)
    
    if( max( [abs(K(i,j) - K(i-1,j)), ...
            abs(K(i,j) - K(i+1,j)), ...
            abs(K(i,j) - K(i,j+1))] )  > 0 )
        borderMat(i,j) = 1;
        im(i,j,:) = [1, zeros(1,size(im,3)-1)];
    end
end

j = size(im,2);
for i = 2:(size(K,1)-1)
    
    if( max( [abs(K(i,j) - K(i-1,j)), ...
            abs(K(i,j) - K(i+1,j)), ...
            abs(K(i,j) - K(i,j-1))] )  > 0 )
        borderMat(i,j) = 1;
        im(i,j,:) = [1, zeros(1,size(im,3)-1)];
    end
end

end