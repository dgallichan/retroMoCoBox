% this script gets run by cluster_submitGRErecon.sh
%
% variable 'rawDataFile' should have been set

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

addpath(getenv('SPM_HOME'));
CLUSTER_LOG_PATH = getenv('CLUSTER_LOG_PATH');
if ~exist('ASPIRE_HOME','var')
    ASPIRE_HOME = [];
end
disp(['ASPIRE_HOME is set as: ' ASPIRE_HOME]);

%%

addpath([retroMoCoPath '/cluster']); % add the cluster subfolder to the path

thisReconPars = getDefaultReconPars('normalRAM');
thisReconPars.rawDataFile = rawDataFile;
thisReconPars.CLUSTER_LOG_PATH = CLUSTER_LOG_PATH;
thisReconPars.ASPIRE_HOME = ASPIRE_HOME;

if exist('bKeepGRAPPArecon','var')
    thisReconPars.bKeepGRAPPArecon = bKeepGRAPPArecon;
end
if exist('outRoot','var')
    thisReconPars.outRoot = outRoot;
end
if exist('parpoolSize','var')
    thisReconPars.parpoolSize = 10;
end



%%

slurmJobID = getenv('SLURM_JOB_ID'); % this will be used for temp folder creation - useful here just to check it works
disp(['SLURM job ID: ' slurmJobID])


reconstructSiemensGREwithFatNavs_cluster(thisReconPars);
