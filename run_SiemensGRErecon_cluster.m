% rawDataFile = ''; % <-- this variable should now be set in calling script for cluster submission

if ~exist(rawDataFile,'file')
    disp(['File ' rawDataFile ' not found, trying locally instead.'])
    rawDataFile = fullfile(script_pwd,rawDataFile);
    if ~exist(rawDataFile,'file')
        disp(['Error: file ' rawDataFile ' not found either'])
    else
        disp(['File ' rawDataFile ' located!'])
    end
end

disp(['Attempting to reconstruct: ' rawDataFile]);

run([getenv('RETROMOCOBOX_HOME') '/addRetroMoCoBoxToPath.m']);
% run([getenv('MIRT_HOME') '/setup.m']); % not needed now that NUFFT from
% MIRT is included in retromocobox directly
addpath(getenv('SPM_HOME'));
CLUSTER_LOG_PATH = getenv('CLUSTER_LOG_PATH');
if ~exist('ASPIRE_HOME','var')
    ASPIRE_HOME = [];
end
disp(['ASPIRE_HOME is set as: ' ASPIRE_HOME]);

%%

addpath([retroMoCoPath '/cluster']); % add the cluster subfolder to the path
%swapDims_xyz = [1 0 1]; % <-- seems to be correct for current FatNav ASPIRE protocol
if ~exist('swapDims_xyz','var')
    swapDims_xyz = [1 0 1]; % <-- seems to be correct for current FatNav ASPIRE protocol
end
if ~exist('bKeepGRAPPArecon','var')
    bKeepGRAPPArecon = 0;
end
if ~exist('outRoot','var')
    outRoot = [];
end
if ~exist('parpoolSize','var')
    parpoolSize = 10;
end
%%

reconstructSiemensGREwithFatNavs_cluster(rawDataFile,'swapDims_xyz',swapDims_xyz,'CLUSTER_LOG_PATH',CLUSTER_LOG_PATH,'bKeepGRAPPArecon',bKeepGRAPPArecon,'outRoot',outRoot,'ASPIRE_HOME',ASPIRE_HOME,'parpoolSize',parpoolSize);
% reconstructSiemensGREwithFatNavs_cluster_debugClusterDependency(rawDataFile,'swapDims_xyz',swapDims_xyz,'CLUSTER_LOG_PATH',CLUSTER_LOG_PATH,'bKeepGRAPPArecon',bKeepGRAPPArecon,'outRoot',outRoot,'ASPIRE_HOME',ASPIRE_HOME,'parpoolSize',parpoolSize);
