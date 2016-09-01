rawDataFile = '/disk2/gallicha_data/temp/meas_MID70_mp2rage_FN600b_FatNav_06mm.dat';

addpath(genpath('~/retroMoCoBox/'))
addpath('~/matlabdownloads/spm12')
run('~/matlabdownloads/fessler/setup.m')

%% if you have plenty of RAM (32 Gb works!):

reconstructSiemensMP2RAGEwithFatNavs(rawDataFile,'bGRAPPAinRAM',1);


%% otherwise do the slower version:
% reconstructMP2RAGEwithFatNavs(rawDataFile,'bGRAPPAinRAM',0);


