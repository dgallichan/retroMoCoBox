function timingReport = reconstructSiemensMP2RAGEwithFatNavs(reconPars)
% function timingReport = reconstructSiemensMP2RAGEwithFatNavs(reconPars)
% 
% Dependencies:
%%%%%%%%%%%%%%%
%
%  %% Required %%:
%
%     SPM 12 (www.fil.ion.ucl.ac.uk/spm) 
%           For coregistration of FatNavs ('spm_realign').
%           This is obviously core to the concept of FatNavs for
%           motion-correction. I have tried using alternative registration
%           software, but SPM appears to have default parameters which work
%           well for FatNavs - being particular sensitive to sub-voxel
%           movements.
%
% Usage:
%%%%%%%%
%
% 'reconPars' - create a struct first using reconPars = getDefaultReconPars()
%               and then modify any parameters as appropriate. 
%
% The only **required** modification to the defaults is to specify the
% rawDataFile itself:
%
%     'reconPars.rawDataFile' - the filename for the raw data you exported from
%                 TWIX (including the full path). Due to the way I chose to
%                 name output files, if you have renamed the datafile it is
%                 required to still contain the text 'MID' followed by a
%                 number and then an underscore ('_') so that multiple
%                 files can be processed from the same directory and not
%                 get overwritten.
%     
%   

%

%  --> overall version notes now pushed to 'retroMocoBoxVersion.m')


%% Validate that the input reconPars contains all required fields

[isValid, errMsg] = validateReconPars(reconPars); 
if ~isValid
    disp(errMsg)
    return
end


%%

reconPars.retroMocoBoxVersion = retroMocoBoxVersion; % update just in case this was old somehow...

rawDataFile = reconPars.rawDataFile;

if ~exist(rawDataFile,'file')
    disp(['Error - raw data file given: "' rawDataFile '" not found']);
    return
end


%%

startTime = clock;

%% Check SPM and Fessler's toolbox are on path

if ~exist('spm.m','file') % could also check version, but that's more effort...
    disp('Error - SPM (ver 12) must be on the path')
    return
end

if ~exist('nufft_init.m','file')
    disp('Error - Fessler toolbox must be on the path for the NUFFT')
    return
end


%%

if reconPars.bFullParforRecon && ~reconPars.bKeepReconInRAM
    disp('Error - you asked for the full parfor option (bFullParforRecon), but not to do the recon in RAM (bKeepReconInRAM)')
    return
end

if reconPars.bFullParforRecon && reconPars.bUseGPU
    disp('Error - you asked for the full parfor option (bFullParforRecon), AND to use GPU (bUseGPU). These are mutually exclusive')
    return
end

%% Create a parpool if it doesn't exist already

tic
isPool = gcp('nocreate');
if isempty(isPool)
    c = parcluster('local');
    if ~isempty(reconPars.parpoolSize)
        c.NumWorkers = reconPars.parpoolSize;
        c.parpool(c.NumWorkers);
    else
        c.parpool();
    end
end
timeToCreateParpool = toc;


%% Disply pie chart of current breakdown of reconstruction time

% Here benchmarked on server allowing parpool size of 64 CPUs and 360 GB
% RAM, testing 600 um and 1000 um data. 25/10/23
% 
% t_TotalTime = timingReport.totalTime; % seconds
% t_ParseRaw = timingReport.parseRawDataFile;
% t_reconFatNavs = timingReport.timingReport_FatNavs.allFatNavs;
% t_SPMrealign = timingReport.timingReport_FatNavs.SPMalignment;
% t_GRAPPArecon = timingReport.timingReport_hostRecon.GRAPPArecon;
% t_NUFFT = timingReport.timingReport_totalTimeApplyMoco;
% t_post = timingReport.postProcessing;
% t_other = t_TotalTime - t_ParseRaw - t_reconFatNavs - t_SPMrealign -t_GRAPPArecon - t_NUFFT - t_post;
% 
% figure(1001)
% % set(gcf,'Position',[    88   415   588   506])
% clf
% pie([t_other t_post t_ParseRaw t_reconFatNavs t_SPMrealign t_GRAPPArecon t_NUFFT ],{'Other','Post-processing','Parse raw data file','Reconstruct FatNavs','SPM realign FatNavs','GRAPPA recon for host','NUFFT'})
% title({'Current breakdown of full reconstruction pipeline on 1 mm data', ['total = ' num2str(t_TotalTime/60,'%.1f') ' mins, running on 64 CPUs with 360 GB RAM'],char(datetime)})
% export_fig('processingTimeBreakdown_1mm_64CPUs_360GB_RAM.png')
% % title({'Current breakdown of full reconstruction pipeline on 600 um data', ['total = ' num2str(t_TotalTime/60,'%.1f') ' mins, running on 64 CPUs with 360 GB RAM'],char(datetime)})
% % export_fig('processingTimeBreakdown_600um_64CPUs_360GB_RAM.png')


          




%% Make raw data object (parses file, but does not load into RAM)

tic
twix_obj = mapVBVD_fatnavs(rawDataFile,'removeOS',1);
timingReport_parseRawDataFile = toc;

if length(twix_obj)>1 % on VE (and presumably VD as well) the raw data typically also has the adjcoilsens data as scan 1, so skip this
    twix_obj = twix_obj{2};
end

if ~isfield(twix_obj,'FatNav')
    disp('Error, no FatNavs found in raw data file!')
    twix_obj % this displays the fields of the twix_obj that are present for comparison/debugging
    return
end

%%


if isempty(reconPars.FatNavRes_mm)
    reconPars.manualFatNavRes = 0;
    nominalB0 = round(twix_obj.hdr.MeasYaps.sProtConsistencyInfo.flNominalB0);
    switch nominalB0
        case 3
            reconPars.FatNavRes_mm = 4;
        case 7
            reconPars.FatNavRes_mm = 2;
        otherwise
            disp(['Error - unexpected field strength of ' nominalB0 'T ...!'])
            return
    end
else
    reconPars.manualFatNavRes = 1;
end

switch reconPars.FatNavRes_mm % in newer version of FatNav sequences, the resolution can be chosen at 2, 4 or 6 mm - with the FOV hard-coded in the sequence itself
    % these choices should also match the corresponding code in processFatNavs_GRAPPA4x4.m
    case {2,4}
        reconPars.FatNav_FOVxyz = [176 256 256]; % FatNav FOV         
        reconPars.FatNav_xyz = reconPars.FatNav_FOVxyz ./ reconPars.FatNavRes_mm;
    case 6
        reconPars.FatNav_FOVxyz = [192 264 264]; % FatNav FOV         
        reconPars.FatNav_xyz = reconPars.FatNav_FOVxyz ./ reconPars.FatNavRes_mm;
        reconPars.FatNav_xyz(3) = 64; % No idea why, but the 6mm data has 64 points in the readout direction for the ACS lines instead of 44...        
end


%% Run the reconstruction on each volume in the raw data

nAve = twix_obj.image.NAve;
nRep = twix_obj.image.NRep;

if nAve > 1

    disp('Multiple Averages detected in raw data file')
    
    for iAve = 1:nAve % can't use parfor here as then it wouldn't work inside each thread
        
        disp(['Processing average ' num2str(iAve)])
        
        thisReconPars = reconPars;
        thisReconPars.iAve = iAve;
        timingReport{iAve} = reconstructSiemensMPRAGEvolume(twix_obj,thisReconPars);
    end
    
else
    if nRep > 1
        disp('Error - handling of data with multiple repetitions not yet implemented...! Please contact gallichand@cardiff.ac.uk for more info')
        return
    else
        timingReport = reconstructSiemensMPRAGEvolume(twix_obj,reconPars);       
    end
end


%%
stopTime = clock;
totalTime = etime(stopTime,startTime)/60/60;
totalTime_hrs = floor(totalTime);
if totalTime_hrs > 0
    totalTime_mins = round(rem(totalTime,totalTime_hrs)*60);
else
    totalTime_mins = round(totalTime*60,2);
end

fprintf('*************************************************************\n')
fprintf('***** reconstructSiemensMP2RAGEwithFatNavs.m completed! *****\n')
fprintf('*************************************************************\n')
fprintf(['Total reconstruction time: ' num2str(totalTime_hrs) ' hours, ' num2str(totalTime_mins) ' mins\n']);

timingReport.totalTime = etime(stopTime,startTime);
timingReport.parseRawDataFile = timingReport_parseRawDataFile;
timingReport.timeToCreateParpool = timeToCreateParpool;
%%


    


end

