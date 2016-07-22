data = ncread('output_huge_test.nc', 'averaged lactate 5 fields');
cylR = ncread('output_huge_test.nc', 'cylR_int');
cylP = ncread('output_huge_test.nc', 'cylP_int');
cylZ = ncread('output_huge_test.nc', 'cylZ_int');

average = squeeze(mean(data, 3));

rCoord = squeeze(cylR(:,:,1));
pCoord = squeeze(cylP(:,:,1));

plotData = squeeze(average(:,:,91));

[x,y] = pol2cart(pCoord,rCoord);

minX = min(min(x));
maxX = max(max(x));

minY = min(min(y));
maxY = max(max(y));

xVec = linspace(minX, maxX,50);
yVec = linspace(minY, maxY,50);

shp = size(plotData);

len = shp(1)*shp(2);

[xMesh, yMesh] = meshgrid(xVec, yVec);
