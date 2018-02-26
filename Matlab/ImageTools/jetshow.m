function [ ] = jetshow( im )
%JETSHOW Summary of this function goes here
%   Detailed explanation goes here

imshow(rescaleIm(im)*64,jet)

end

