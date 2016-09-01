function sn(data,filename,voxdims, origin)
% function sn(data,filename,voxdims,origin)
%
% Save a nifti file with an empty header

if nargin < 3
    voxdims = [1 1 1];
end

if nargin < 4
    szData = size(data);
    origin = szData(1:3)/2+1;
end

temphdr = make_nii(data,voxdims,origin);

if ~strcmp(filename(end-3:end),'.nii')
    filename = [filename '.nii'];
end
save_nii(temphdr,filename);
