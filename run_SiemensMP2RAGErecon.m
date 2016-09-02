% Obviously this needs to point to your data (Siemens VB17 MP2RAGE with 3D-FatNavs)
rawDataFile = '/home/gallicha/temp/autoTransferTemp/tempRaw/meas_MID36_mp2rage_FN600b_FatNav_06mm.dat';

% And wherever you installed the Retro-MoCo-Box
addpath(genpath('~/retroMoCoBox/'))

% And SPM 12
addpath('~/matlabdownloads/spm12')

% And Prof. J. Fessler's Matlab toolbox (http://web.eecs.umich.edu/~fessler/code/) for the NUFFT
run('~/matlabdownloads/fessler/setup.m')

%% Fastest option: if you have loads of RAM (tested with 12 CPUs and 96 Gb of RAM on both 1 mm and 600 um data)

reconstructSiemensMP2RAGEwithFatNavs(rawDataFile,'bGRAPPAinRAM',1,'bKeepReconInRAM',1,'bFullParforRecon',1);

%% if you have plenty of RAM:

% reconstructSiemensMP2RAGEwithFatNavs(rawDataFile,'bGRAPPAinRAM',1);


%% otherwise do the slower version:
% reconstructMP2RAGEwithFatNavs(rawDataFile,'bGRAPPAinRAM',0);


