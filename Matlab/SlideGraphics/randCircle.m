function [ circle ] = randCircle( center, radius, numPoints )
%RANDCIRCLE Summary of this function goes here
%   Detailed explanation goes here

circle = (rand(numPoints,2) - [0.5 0.5])*2;

rownorms = sum(circle.^2,2);
circle(rownorms > 1,:) = [];

circle = circle*radius + center;

end

