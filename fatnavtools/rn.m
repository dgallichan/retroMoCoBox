function out = rn(filename)
% function out = rn(filename)
%
% Read in a nifti file without 'touching' it...

inName = filename;
if ~exist(filename,'file')
    if strcmp(filename(end-2:end),'nii')        
        filename = [filename '.gz'];
    else
        filename = [filename '.nii'];
        if ~exist(filename,'file')
            filename = [filename '.gz'];
        end
    end
end

if ~exist(filename,'file')
    disp(['Error, cannot find file: ' inName ])
    return;
end
    
temp = load_untouch_nii(filename);
out = double(temp.img);