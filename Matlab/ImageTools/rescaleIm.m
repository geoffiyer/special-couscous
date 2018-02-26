function [ im ] = rescaleIm( im )
%RESCALEIM Summary of this function goes here
%   Detailed explanation goes here

if(size(im,3) > 1)
   themax = max(max(max(im))); 
   themin = min(min(min(im)));
   im = (im - themin)/(themax - themin);
else
   themax = max(max(im)); 
   themin = min(min(im));
   im = (im - themin)/(themax - themin);
end

end

