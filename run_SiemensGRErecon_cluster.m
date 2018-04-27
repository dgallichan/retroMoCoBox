% rawDataFile = ''; % <-- this variable should now be set in calling script for cluster submission

disp(['Attempting to reconstruct: ' rawDataFile]);

run([getenv('RETROMOCOBOX_HOME') '/addRetroMoCoBoxToPath.m']);
run([getenv('MIRT_HOME') '/setup.m']);
addpath(getenv('SPM_HOME'));
CLUSTER_LOG_PATH = getenv('CLUSTER_LOG_PATH');

%%

addpath([moCoBoxPath '/cluster']); % add the cluster subfolder to the path
swapDims_xyz = [1 0 1]; % <-- seems to be correct for current FatNav ASPIRE protocol

%%

reconstructSiemensGREwithFatNavs_cluster(rawDataFile,'swapDims_xyz',swapDims_xyz,'CLUSTER_LOG_PATH',CLUSTER_LOG_PATH);
