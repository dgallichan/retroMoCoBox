% add retroMoCoBox subfolders to the path (without including all the git
% subfolders!)

retroMoCoPath = which('addRetroMoCoBoxToPath.m');
[retroMoCoPath, ~] = fileparts(retroMoCoPath);

addpath(retroMoCoPath);
addpath([retroMoCoPath '/export_fig']);
addpath([retroMoCoPath '/fatnavtools']);
addpath([retroMoCoPath '/generaltools']);
addpath([retroMoCoPath '/images']);
addpath([retroMoCoPath '/niftitools']);
addpath([retroMoCoPath '/mirt_nufft']);
addpath([retroMoCoPath '/mirt_nufft/mex']);
addpath([retroMoCoPath '/mirt_nufft/utilities']);
if exist([retroMoCoPath '/mapVBVD_20160905'],'dir')
    addpath([retroMoCoPath '/mapVBVD_20160905']);
end