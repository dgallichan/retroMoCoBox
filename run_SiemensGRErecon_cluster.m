% rawDataFile = ''; % <-- this variable should now be set in calling script for cluster submission

disp(['Attempting to reconstruct: ' rawDataFile]);

run([getenv('RETROMOCOBOX_HOME') '/addRetroMoCoBoxToPath.m']);
run([getenv('MIRT_HOME') '/setup.m']);
addpath(getenv('SPM_HOME'));
CLUSTER_LOG_PATH = getenv('CLUSTER_LOG_PATH');

%%

addpath([retroMoCoPath '/cluster']); % add the cluster subfolder to the path
%swapDims_xyz = [1 0 1]; % <-- seems to be correct for current FatNav ASPIRE protocol
if ~exist('swapDims_xyz','var')
    swapDims_xyz = [1 0 1]; % <-- seems to be correct for current FatNav ASPIRE protocol
end
if ~exist('bKeepGRAPPArecon','var')
    bKeepGRAPPArecon = 0;
end

%%

reconstructSiemensGREwithFatNavs_cluster(rawDataFile,'swapDims_xyz',swapDims_xyz,'CLUSTER_LOG_PATH',CLUSTER_LOG_PATH,'bKeepGRAPPArecon',bKeepGRAPPArecon);

