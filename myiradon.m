function BPI = myiradon(sinogram,thetas)

numOfParallelProjections = size(sinogram,1);
numOfAngularProjections  = length(thetas); 
thetas = (pi/180)*thetas;

BPI = zeros(numOfParallelProjections,numOfParallelProjections);
midindex = floor(numOfParallelProjections/2) + 1;

[xCoords,yCoords] = meshgrid(ceil(-numOfParallelProjections/2):ceil(numOfParallelProjections/2-1));

for i = 1:numOfAngularProjections
    rotCoords = round(midindex + xCoords*sin(thetas(i)) + yCoords*cos(thetas(i)));
    indices   = find((rotCoords > 0) & (rotCoords <= numOfParallelProjections));
    newCoords = rotCoords(indices);    
    BPI(indices) = BPI(indices) + sinogram(newCoords,i)./numOfAngularProjections;
    drawnow
end