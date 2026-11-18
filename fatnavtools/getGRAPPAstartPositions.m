function [startX, startY] = getGRAPPAstartPositions(data2D,grapKernel)
% function [startX, startY] = getGRAPPAstartPositions(data2D,grapKernel)
%
% based on the data and the kernel, find where to start applying the kernel

if ndims(data2D)>2
    error('data2D must be 2D!')
end

dataMask = zeros(size(data2D));
dataMask(data2D~=0) = 1;

[gx,gy] = size(grapKernel);
[sx,sy] = size(data2D);

shifts_x = 0:gx-1;
shifts_y = 0:gy-1;

shiftScores = zeros(gx-1,gy-1);

thisGrapKernel = grapKernel;
thisGrapKernel(thisGrapKernel==0.5) = 0; % remove the data to fill in

kernelMask = repmat(thisGrapKernel(1:end-1,1:end-1),ceil(sx/(gx-1)),ceil(sy/(gy-1)));
kernelMask = kernelMask(1:sx,1:sy); % cut to size

for iX = 1:gx-1
    for iY = 1:gy-1
        thisKernelMask = circshift(kernelMask,[shifts_x(iX) shifts_y(iY)]);
        shiftScores(iX,iY) = sum(reshape(thisKernelMask.*dataMask,sx*sy,1));
    end
end

% get the indices of the max closest to 0,0
[iX,iY] = find(shiftScores==max(shiftScores(:)));
shiftMags = sqrt(shifts_x(iX).^2 + shifts_y(iY).^2);
iMinShiftMag = find(shiftMags==min(shiftMags),1);

startX = iX(iMinShiftMag);
startY = iY(iMinShiftMag);
