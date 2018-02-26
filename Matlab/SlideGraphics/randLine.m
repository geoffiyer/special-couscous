function [ line ] = randLine( startpoint, endpoint, width, numpoints )
% input: startpoint - 1x2 vector with x,y coords of start
%        endpoint   - 1x2 vector with x,y coords of end
%        width      - scalar with max width
%        numpoints  - scalar with number of points
%output: line       - numpointsx2 matrix. x,y coords of each point.

lineVec = endpoint-startpoint;
perpVec = [-lineVec(2) lineVec(1)];
lineLength = norm(lineVec);

line = rand(numpoints,1)*lineVec + startpoint;
variance = (rand(numpoints,1)-0.5)*perpVec*width/lineLength;

line = line + variance;

end

