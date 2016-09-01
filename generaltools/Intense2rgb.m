function out = Intense2rgb(data,cmap,clims)
% function out = Intense2rgb(data,cmap,clims)
% take 2d data and return rgb image, based on given colormap (gray is
% default)

if nargin < 3
    clims = [min(data(:)) max(data(:))];
end

if ( nargin < 2 | isempty(cmap) )
    cmap = repmat(linspace(0,1,64)',1,3);
end

cmapLength = size(cmap,1);

data(data>clims(2)) = clims(2);
data(data<clims(1)) = clims(1);

data = (data-clims(1))./(clims(2)-clims(1));

mask = ones(size(data));

dataMat = vols2matrix(data,mask);

outMat = zeros(size(dataMat,1),3);
outMat(:,1) = interp1(linspace(0,1,cmapLength),cmap(:,1),dataMat);
outMat(:,2) = interp1(linspace(0,1,cmapLength),cmap(:,2),dataMat);
outMat(:,3) = interp1(linspace(0,1,cmapLength),cmap(:,3),dataMat);

out = reshape(matrix2vols(outMat,mask),[size(data,1) size(data,2) 3]);
