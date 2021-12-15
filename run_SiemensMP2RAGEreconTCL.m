% Obviously this needs to point to your data (Siemens MP(2)RAGE with 3D-FatNavs)
rawDataFile = '/Volumes/WD2TB_Data2/data/2018_05_10_TracInnovTest/meas_MID00045_FID34860_MPRAGE_FatNavs_1mm_iso_movement.dat';

% And this needs to point to the folder containing the TCL data from the
% same session:
TCLdir = '/Volumes/WD2TB_Data2/data/2018_05_10_TracInnovTest/DataTCL/2018-05-10_MPRAGE';

% And the TCL Suite needs to be on the path too
addpath ~/matlab/matlabdownloads/tracSuite/

% And wherever you installed the Retro-MoCo-Box
run('~/Documents/code/retroMoCoBox/addRetroMoCoBoxToPath.m')

% And the Michigan Image Reconstruction Toolbox (MIRT) (http://web.eecs.umich.edu/~fessler/code/) for the NUFFT
% run('~/matlab/matlabdownloads/mirt/setup.m') % mirt_nufft now added to retromocobox

TCLtimeOffset_ms = -82229 + 900; % -82229 is from manual logging that day, 900 is observed additional offset in other data

%% Fastest option: if you have loads of RAM (tested with 12 CPUs and 96 Gb of RAM on both 1 mm and 600 um data)

reconstructSiemensMP2RAGEwithTCL(rawDataFile,TCLdir,'TCLtimeOffset_ms',TCLtimeOffset_ms,'bGRAPPAinRAM',1,'bKeepReconInRAM',1,'bFullParforRecon',1);

%% if you have plenty of RAM:

% reconstructSiemensMP2RAGEwithTCL(rawDataFile,TCLdir,'bGRAPPAinRAM',1);


%% otherwise do the slower version:

% reconstructSiemensMP2RAGEwithTCL(rawDataFile,TCLdir,'bGRAPPAinRAM',0);


