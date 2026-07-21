clear
close all


% Obviously this needs to point to your data (Siemens MP(2)RAGE with 3D-FatNavs)
rawDataFile = '/home/scedg10/mydata/testdata/meas_MID132_mp2rage_FN600b_FatNav_06mm.dat';

% And wherever you installed the Retro-MoCo-Box
run('~/retroMoCoBox/addRetroMoCoBoxToPath.m')

% And SPM 12
addpath('/cubric/software/spm.versions/spm12')

% You don't have to specify the following (default if empty is to use the
% same folder where the .dat is located) - a new folder based on the MID
% number of the .dat file is created here
outRootFolder = '/home/scedg10/mydata/testdata/testrecons';

%% Fastest option: if you have loads of RAM (tested with 12 CPUs and 96 Gb of RAM on both 1 mm and 600 um data)
%
% (Note that there is a bit of a 'trick' available - if you have a
% high-power server with many CPUs and high RAM, but still a very large
% dataset that won't fit into RAM enough times over for full parfor on all
% CPUs, you can set the parpool size down to say ~8-10 and it might still
% work reasonably fast then! - use e.g. theseReconPars.parpoolSize = 8;

thisReconPars = getDefaultReconPars('veryhighRAM');
thisReconPars.rawDataFile = rawDataFile; % required

thisReconPars.outRoot = outRootFolder; % optional
thisReconPars.outFolderPrefix = 'veryhighRAM'; % optional - allows multiple recons to exist with same MID in same folder

reconstructSiemensMP2RAGEwithFatNavs(thisReconPars);


%% if you have plenty of RAM:

thisReconPars = getDefaultReconPars('highRAM');
thisReconPars.rawDataFile = rawDataFile;

thisReconPars.outRoot = outRootFolder; % optional
thisReconPars.outFolderPrefix = 'highRAM'; % optional - allows multiple recons to exist with same MID in same folder

reconstructSiemensMP2RAGEwithFatNavs(thisReconPars);

%% otherwise do the slower version:

thisReconPars = getDefaultReconPars('normalRAM');
thisReconPars.rawDataFile = rawDataFile;

thisReconPars.outRoot = outRootFolder; % optional
thisReconPars.outFolderPrefix = 'normalRAM'; % optional - allows multiple recons to exist with same MID in same folder

reconstructSiemensMP2RAGEwithFatNavs(thisReconPars);


