% rawDataFile = '';  % <-- this must already be set before calling this script

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
%run([getenv('MIRT_HOME') '/setup.m']);
addpath(getenv('SPM_HOME'));

%%
dir(rawDataFile)
%%

thisReconPars = getDefaultReconPars('veryhighRAM');
thisReconPars.rawDataFile = rawDataFile;
thisReconPars.outRoot = outRoot;
thisReconPars.FatNavRes_mm = FatNavRes_mm;

reconstructSiemensMP2RAGEwithFatNavs(thisReconPars);


