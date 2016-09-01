function [dat2]=vols2matrix(data,mask)
% [dat2,ind]=vols2matrix(data,mask);
%
% takes a 3d/4d volume and mask (2d/3d) and returns the 2d matrix 
% (space x time) 

maskdims = size(mask);
mask=reshape(mask,numel(mask),1)'>0;
dims = size(data);
if length(dims) == length(maskdims) % i.e. time dimension == 1
    dims(end+1) = 1;
end
dat2=reshape(data,numel(mask),dims(end))';
dat2=dat2(:,mask)';
