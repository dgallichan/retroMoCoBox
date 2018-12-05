function reconstructSiemensMP2RAGEwithFatNavs(rawDataFile,varargin)
% function reconstructSiemensMP2RAGEwithFatNavs(rawDataFile,varargin)
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
%     Michigan Image Reconstruction Toolbox (MIRT)
%     (http://web.eecs.umich.edu/~fessler/code/index.html)
%     from the group of Prof. J Fessler
%           This is used to perform the NUFFT 3D gridding operation used to deal
%           with the non-Cartesian k-space sampling after rotations have
%           been applied.
%
% Usage:
%%%%%%%%
%
% 'rawDataFile' - this is the filename for the raw data you exported from
%                 TWIX (including the full path). Due to the way I chose to
%                 name output files, if you have renamed the datafile it is
%                 required to still contain the text 'MID' followed by a
%                 number and then an underscore ('_') so that multiple
%                 files can be processed from the same directory and not
%                 get overwritten.
%
%  Optional arguments:-
%      'outRoot' - this gets used as the output folder for all processing
%                 (and the default root for the temporary folder which may
%                 generate several Gb of temporary files along the way - 
%                 - override this with 'tempRoot' option). The default 
%                 'outRoot' is to use the same folder as the raw data file.
%    
%      'tempRoot' - location for creating temporary files - can be many Gb.
%                   Default is to use the same location as 'outRoot'.
%
%      'bLinParSwap' - set this to '1' to indicate that the 'LIN/PAR swap'
%                     option was chosen in the MP2RAGE sequence. This
%                     alters which direction in k-space the FatNavs
%                     correspond to.
%
%      'bGRAPPAinRAM' - set this to '1' to perform the GRAPPA recon (if
%                       necessary!) for the host sequence entirely in RAM. This
%                       can be quite a lot faster - but requires sufficient
%                       RAM...! 
%                       (Note that the temporary variables for applying the
%                       MoCo are currently done outside of RAM as this
%                       is currently not the bottleneck for that part of
%                       the recon)
%
%      'bKeepGRAPPArecon' - set this to '1' to prevent deleting of the
%                           GRAPPA recon of the host data (which would be
%                           deleted by default). 
%
%      'bKeepReconInRAM' - set this to '1' to keep all the reconstructed
%                          files (INV1, INV2 and UNI, all with and without
%                          correction) in RAM, or default is to use a
%                          temporary file.
%
%      'bFullParforRecon' - set this to '1' to enable parfor over the NUFFT
%                           loop. For this to work, bKeepReconInRAM must
%                           also be set - and depending how many CPUs you
%                           have available on your matlabpool, you may need
%                           to have LOADS of RAM...  But it is much faster!
%       
%      'coilCombineMethod' - In MP2RAGE the default method is to match what
%                            I believe is done on the scanner - to weight
%                            each separately calculated UNI image by the
%                            square of the INV2 image for that coil.
%                            Slightly better results might be obtained for
%                            low SNR data by using a lower-res version of
%                            the INV2 image. Set this option to 'lowres' to
%                            try this. 
%                            For INV1 and INV2 the default is to combine
%                            the images using root sum-of-squares
%                            - this is not optimal, so the 'lowres' will
%                            also try to use a low-res version of INV2 to
%                            combine both. This corresponds to the coil
%                            combination method of Bydder et al, MRM 2002,
%                            47:539-548 - but as there is a smoothness
%                            parameter which needs tuning - and that it can
%                            also lead to signal voids, I haven't yet felt
%                            confident enough to make this the default
%                            processing.
%
%       'FatNavRes_mm' - the spatial resolution of the acquired FatNavs, 
%                        specified in mm. The current implementation of the
%                        pulse sequence does not allow this information to
%                        be stored in the raw data, so this must be entered
%                        here manually - or use the default of 2 mm.
%                        Whatever you specify here, the code does still
%                        assume that the acceleration for the FatNavs was
%                        4x4 and that the FOV was 176x256x256mm.
%
%       'swapDims_xyz' - 3-component row-vector to determine whether to
%                        reverse each of x, y and z directions. Default is
%                        [0 0 1] which seems to work for a lot of parameter
%                        sets, but not all...  Check the orientation
%                        checker feature in the HTML!
%
%       'bZipNIFTIs' - Use '1' to apply gzip at the end to all the NIFTI
%                     files (default), otherwise just leave them uncompressed.
%
%       'bKeepFatNavs' - Use '1' to keep all the reconstructed FatNavs,
%                       otherwise delete that folder when finished
%                       (default).
%
%       'bKeepPatientInfo' - Use '1' to keep sensitive info from raw data header  
%                           in HTML output (default). Use '0' to anomyize completely
%                           and use a string e.g. '0019' to insert the ID from
%                           another database.
%
%
%
%     
%   
% Matlab tools which are included (with 'assumed' permission, as I collected them online):
%
%     mapVBVD (from Philip Ehses) 
% *** Note that this code is presumably 'Siemens-sensitive' as it is    ***
% *** not freely available online, but only from the Siemens user forum.***
% *** Consequently the FatNavs recon code must also be considered       ***
% *** 'Siemens-sensitive' while it contains this code                   ***
% *** I have not included this part in the Github repository - please   ***
% *** email me if you would like it.                                    ***
%           For reading in the Siemens raw data format - and able to handle
%           very large datasets elegantly. 
%           I made small changes to the code to allow handling of the
%           FatNav data.
%
%     NIFTI Tools (from Jimmy Shen)
%           For reading in and saving out in the .nii NIFTI format. 
%           This code is provided here unaltered.            
%
%     GCC - coil compression (from Miki Lustig)
%           We work mostly with a 32-channel head coil, and so massive
%           speedups are possible when using coil compression. This code
%           from Miki Lustig's website implements the method described in
%           Zhang et. al MRM 2013;69(2):571-82.
%           For this code I have directly taken snippets and inserted them
%           into performHostGRAPPArecon.m
% 
%     export_fig (from Oliver Woodford)
%           Useful for making figures that you save from Matlab look nice!
%           :)
%
%     process_options (from Mark A. Paskin)
%           Useful for sending optional inputs to this master file.         
%
%     subplot1 (from Eran O. Ofek)
%           Nicer than Matlab's own way of doing subplots
%
% 
% Optional:
%     ImageMagick (www.imagemagick.org)
%           For making animated GIFs of results ('convert') in the html
%           report
%
%     FSL (www.fmrib.ox.ac.uk/fsl) - tested with v5.0
%           Used for brain extraction ('bet') in order to make MIP views 
%           of INV2 images (which show bright arteries) more interpretable
%
%
%
%
%
% To do:
%%%%%%%%
%
%  Coil compression for the FatNavs
%              - this is now implemented, but preliminary tests have shown
%                it doesn't work well... One problem is the weights
%                themselves for the sparse image - the other problem is
%                that acceleration of 16 is quite a lot...
%
%  Coil compression when no GRAPPA used
%              - this might make sense to still be able to speed up the
%                application of the retrospective motion-correction, but as
%                most people are probably scanning with GRAPPA in the host
%                sequence anyway, I haven't implemented this yet.
%
%
% -------------------------------------------------------------------------
% reconstructSiemensMP2RAGEwithFatNavs.m, 
%   v0.1 - daniel.gallichan@epfl.ch - June 2015
%   v0.2 -    -- February 2016 --
%        - added option to specify resolution of FatNavs 
%        - also added option to specify swapDims_xyz 
%        - added Patient info to HTML output
%        - added option to zip the NIFTIs at the end
%        - added option to keep FatNavs (and changed default behaviour to delete them)
%   v0.3 -   -- April 2016 -- trying to speed things up
%        - changed default oversampling from 2 down to 1.375 for NUFFT
%        - use parfor in NUFFT
%        - use parfor in SPM registration
%        - switch to using SPM12 (previously SPM8)
%        - now have option to do the MP2RAGE recon combination entirely in RAM
%   v0.4 -   -- August 2016
%        - changed default NUFFT oversampling from 1.375 to 1.5 (to reduce aliasing artifact)
%   v0.5 -   -- September 2016
%        - updated to latest version of Philipp Ehses' mapVBVD software (from 18/9/15)
%        - now renamed things to make it part of the 'RetroMoCoBox' and put
%          on Github
%        - Added the bFullParforRecon option to really speed things up if
%          you have enough CPUs and RAM available
%        - changed name of bKeepRecoInRAM to bKeepReconInRAM for consistency
%        - output files no longer start with 'a_host_'
%        - fit to versioning for the whole of 'RetroMoCoBox' as 0.5.0
%        
%   0.5.1 -  -- September 2016
%         - Fixed bug in handling of data acquired without GRAPPA
%         - Fixed bug in handling of data acquired with different orientations 
%
%   0.6.0 -  -- February 2017 - new contact email: gallichand@cardiff.ac.uk
%         - *CHANGED* handling of motion estimates - now average temporal
%           neighbours
%         - Add support for VD/VE data
%         - *RENAMED* output 'uniImage' and 'uniImage_corrected' to 'UNI'
%           and 'UNI_corrected' as seems to match INV1 and INV2 better
%         - Add animated GIF with zoom of front of brain where changes are
%           likely to be most noticeable and put in HTML
%
%   0.6.1 - -- May 2017
%         - Added Torben's option to anonymize data
%
%   0.6.2 -   -- July 2017
%         - Improved 'parfor' handling of data without second inversion
%           time
%         - Added manual discard of channels for HeadNeck_64 coil which
%           have a lot of signal in the neck
%
%   0.6.3 - -- August 2017
%         - Include the FatNav resolution in the HTML (and display to
%           screen)
%         - handle the case where 'PatientName' becomes 'tPatientName' for
%           no apparent reason
%
%   0.6.4 - -- Sep 2017
%         - Automatically set FatNavRes_mm based on field strength (7T - 2mm, 3T - 4mm)         
%
%   0.7.0 - Feb 2018
%         - Create sub-function reconstructSiemensVolume.m so that multiple
%           averages or repetitions can be handled directly in this code.
%           NB. Handling this properly would involve updated the sequence
%           code run on the scanner because it currently doesn't label the
%           FatNavs properly beyond the first volume of the host. This
%           feature probably won't be used enough to make that worthwhile
%           though...
%         - WARNING - I took this opportunity to cleanup the names of some
%           of the input options (always putting a 'b' for boolean in front
%           of logical options) so please check your calling code!
%
%
%   0.7.1 - -- Nov 2018
%         - Added new input flag 'bKeepComplexImageData' to allow saving
%           out of the complex data per coil as MATLAB files - before and
%           after application of MoCo.

reconPars.retroMocoBoxVersion = '0.7.1dev'; % put this into the HTML for reference
reconPars.rawDataFile = rawDataFile;

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

[reconPars.outRoot, reconPars.tempRoot, reconPars.bLinParSwap, reconPars.bGRAPPAinRAM, reconPars.bKeepGRAPPArecon, reconPars.bKeepReconInRAM, reconPars.bFullParforRecon,...
    reconPars.coilCombineMethod, reconPars.FatNavRes_mm, reconPars.swapDims_xyz, reconPars.bZipNIFTIs, reconPars.bKeepFatNavs,reconPars.bKeepPatientInfo,...
    reconPars.bKeepComplexImageData] = process_options(varargin,...
    'outRoot',[],'tempRoot',[],'bLinParSwap',0,'bGRAPPAinRAM',0,'bKeepGRAPPArecon',0,'bKeepReconInRAM',0,...
    'bFullParforRecon',0,'coilCombineMethod','default','FatNavRes_mm',[],'swapDims_xyz',[0 0 1],'bZipNIFTIs',1,'bKeepFatNavs',0,'bKeepPatientInfo',1,...
    'bKeepComplexImageData',0);


%%

if reconPars.bFullParforRecon && ~reconPars.bKeepReconInRAM
    disp('Error - you asked for the full parfor option (bFullParforRecon), but not to do the recon in RAM (bKeepReconInRAM)')
    return
end

%% Disply pie chart of current breakdown of reconstruction time

% % Here benchmarked on server allowing parpool size of 12 CPUs and 96 Gb
% % RAM, testing 600 um data
% 
% t_TotalTime = 24;
% t_ParseRaw = 86/60;
% t_reconFatNavs = 332/60;
% t_SPMrealign = 61/60;
% t_GRAPPArecon = 6;
% t_NUFFT = 8;
% t_other = t_TotalTime - t_ParseRaw - t_reconFatNavs - t_SPMrealign - t_NUFFT;
% 
% figure(1001)
% set(gcf,'Position',[    88   415   588   506])
% clf
% pie([t_other t_ParseRaw t_reconFatNavs t_SPMrealign t_GRAPPArecon t_NUFFT],{'Other','Parse raw data file','Reconstruct FatNavs','SPM realign FatNavs','GRAPPA recon for host','NUFFT'})
% title({'Current breakdown of full reconstruction pipeline on 600 um data', ['total = ' num2str(t_TotalTime) ' mins, running on 12 CPUs with 96 Gb RAM'],char(datetime)})
% fontScale(1.4)
% export_fig('processingTimeBreakdown.png')


          

%%

startTime = clock;


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
        timingReport{iAve} = reconstructSiemensVolume(twix_obj,thisReconPars);
    end
    
else
    if nRep > 1
        disp('Error - handling of data with multiple repetitions not yet implemented...! Please contact gallichand@cardiff.ac.uk for more info')
        return
    else
        timingReport = reconstructSiemensVolume(twix_obj,reconPars);       
    end
end


%%
stopTime = clock;
totalTime = etime(stopTime,startTime)/60/60;
totalTime_hrs = floor(totalTime);
if totalTime_hrs > 0
    totalTime_mins = round(rem(totalTime,totalTime_hrs)*60);
else
    totalTime_mins = round(totalTime*60);
end

fprintf('*************************************************************\n')
fprintf('***** reconstructSiemensMP2RAGEwithFatNavs.m completed! *****\n')
fprintf('*************************************************************\n')
fprintf(['Total reconstruction time: ' num2str(totalTime_hrs) ' hours, ' num2str(totalTime_mins) ' mins\n']);

%%


    


end

