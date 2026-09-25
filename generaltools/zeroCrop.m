function [out, lims] = zeroCrop(data)
% up to 3D, crop out edge zeros

lims.xMin = find(any(any(data(:,:,:,1),2),3),1,'first');
lims.xMax = find(any(any(data(:,:,:,1),2),3),1,'last');

lims.yMin = find(squeeze(any(any(data(:,:,:,1),1),3)),1,'first');
lims.yMax = find(squeeze(any(any(data(:,:,:,1),1),3)),1,'last');

lims.zMin = find(squeeze(any(any(data(:,:,:,1),1),2)),1,'first');
lims.zMax = find(squeeze(any(any(data(:,:,:,1),1),2)),1,'last');

szData = size(data);
out = data(lims.xMin:lims.xMax,lims.yMin:lims.yMax,lims.zMin:lims.zMax,:);
szOut = size(out);
out = reshape(out,[szOut(1:3) szData(4:end)]);

