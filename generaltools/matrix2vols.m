function out = matrix2vols(data,mask)
%function out = matrix2vols(data,mask)
%
% takes 2d matrix [space x time] and 2d/3d mask [x y (z)]
% and returns 3d/4d data [x y (z) time]

mask2= reshape(mask,numel(mask),1)>0;
out=zeros(numel(mask),size(data,2));
out(mask2,:)=data;

out=reshape(out,[size(mask) size(data,2)]);
