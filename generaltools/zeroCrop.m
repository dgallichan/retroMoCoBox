function out = zeroCrop(data)
% up to 3D, crop out edge zeros

xMin = find(any(any(data(:,:,:,1),2),3),1,'first');
xMax = find(any(any(data(:,:,:,1),2),3),1,'last');

yMin = find(squeeze(any(any(data(:,:,:,1),1),3)),1,'first');
yMax = find(squeeze(any(any(data(:,:,:,1),1),3)),1,'last');

zMin = find(squeeze(any(any(data(:,:,:,1),1),2)),1,'first');
zMax = find(squeeze(any(any(data(:,:,:,1),1),2)),1,'last');

szData = size(data);
out = data(xMin:xMax,yMin:yMax,zMin:zMax,:);
szOut = size(out);
out = reshape(out,[szOut(1:3) szData(4:end)]);

