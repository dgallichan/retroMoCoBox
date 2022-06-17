% rawDataFile = '';  % <-- this must already be set before calling this script

disp(['Attempting to reconstruct: ' rawDataFile]);

run([getenv('RETROMOCOBOX_HOME') '/addRetroMoCoBoxToPath.m']);
%run([getenv('MIRT_HOME') '/setup.m']);
addpath(getenv('SPM_HOME'));

%% Haven't yet worked out a good way to feed these options from the command line launch for the cluster
if ~exist('FatNavRes_mm','var')
    FatNavRes_mm = 2;
end
if ~exist('swapDims_xyz','var')
    swapDims_xyz = [0 0 1]; % seems to be most common...
end
if ~exist('bKeepGRAPPArecon','var')
    bKeepGRAPPArecon = 0;
end
if ~exist('bLinParSwap','var')
    bLinParSwap = 0;
end
if ~exist('bFullParforRecon','var')
    bFullParforRecon = 1;
end
if ~exist('bKeepComplexImageData','var')
    bKeepComplexImageData = 0;
end
if ~exist('outRoot','var')
    outRoot = [];
end

%%
dir(rawDataFile)
%%

reconstructSiemensMP2RAGEwithFatNavs(rawDataFile,'FatNavRes_mm',FatNavRes_mm,...
    'bGRAPPAinRAM',1,'bKeepReconInRAM',1,'bFullParforRecon',bFullParforRecon,'swapDims_xyz',swapDims_xyz,'bKeepFatNavs',0,...
    'bKeepGRAPPArecon',bKeepGRAPPArecon,'bLinParSwap',bLinParSwap,'bKeepComplexImageData',bKeepComplexImageData,...
    'outRoot',outRoot);


