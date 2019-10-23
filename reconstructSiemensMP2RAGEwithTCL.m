function reconstructSiemensMP2RAGEwithTCL(rawDataFile,TCLdir,varargin)
% function reconstructSiemensMP2RAGEwithTCL(rawDataFile,TCLdir,varargin)
% 
% Dependencies:
%%%%%%%%%%%%%%%
%
%  %% Required %%:
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
% 'TCLdir' - this is the folder containing the TCL (TracInnovations
%            Tracoline system) logging information from the same scan
%            session.
%
%
%  Optional arguments:-
%      'TCLtimeOffset_ms' - manually specify time offset between TCL clock
%                           and raw data of MP2RAGE scan
%
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
%       'swapDims_xyz' - 3-component row-vector to determine whether to
%                        reverse each of x, y and z directions. Default is
%                        [0 0 1] which seems to work for a lot of parameter
%                        sets, but not all...  Check the orientation
%                        checker feature in the HTML!
%
%       'bZipNIFTIs' - Use '1' to apply gzip at the end to all the NIFTI
%                     files (default), otherwise just leave them uncompressed.
%
%       'bKeepPatientInfo' - Use '1' to keep sensitive info from raw data header  
%                           in HTML output (default). Use '0' to anomyize completely
%                           and use a string e.g. '0019' to insert the ID from
%                           another database.
%
%       'GRAPPAlambda' - regularization for GRAPPA weights for MPRAGE
%                        recon. Default is zero, but for large motion may
%                        want to use 1e-2 or so.
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
%           FatNav data. For the TCL data, certain properties of my modified 
%           version are also used, despite having no FatNavs.
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
% -- July 2018 - gallichand@cardiff.ac.uk
% -- modification of original reconstructSiemensMP2RAGEwithFatNavs.m to be
%    adapted for motion paramters from the Tracoline TCL system

reconPars.retroMocoBoxVersion = retroMocoBoxVersion; % put this into the HTML for reference
reconPars.rawDataFile = rawDataFile;
reconPars.TCLdir = TCLdir;

%% Check MIRT toolbox is on path

if ~exist('nufft_init.m','file')
    disp('Error - MIRT toolbox must be on the path for the NUFFT')
    return
end

%%

[reconPars.TCLtimeOffset_ms, reconPars.outRoot, reconPars.tempRoot, reconPars.bLinParSwap, reconPars.bGRAPPAinRAM, reconPars.bKeepGRAPPArecon, reconPars.bKeepReconInRAM, reconPars.bFullParforRecon,...
    reconPars.coilCombineMethod, reconPars.swapDims_xyz, reconPars.bZipNIFTIs, reconPars.bKeepPatientInfo, reconPars.GRAPPAlambda] = process_options(varargin,...
    'TCLtimeOffset_ms',0,'outRoot',[],'tempRoot',[],'bLinParSwap',0,'bGRAPPAinRAM',0,'bKeepGRAPPArecon',0,'bKeepReconInRAM',0,...
    'bFullParforRecon',0,'coilCombineMethod','default','swapDims_xyz',[0 0 1],'bZipNIFTIs',1,'bKeepPatientInfo',1,'GRAPPAlambda',0);


%%

if reconPars.bFullParforRecon && ~reconPars.bKeepReconInRAM
    disp('Error - you asked for the full parfor option (bFullParforRecon), but not to do the recon in RAM (bKeepReconInRAM)')
    return
end
          

%%

startTime = clock;


%% Make raw data object (parses file, but does not load into RAM)

tic
twix_obj = mapVBVD_fatnavs(rawDataFile,'removeOS',1);
timingReport_parseRawDataFile = toc;

if length(twix_obj)>1 % on VE (and presumably VD as well) the raw data typically also has the adjcoilsens data as scan 1, so skip this
    prescan_twix_obj = twix_obj{1};
    twix_obj = twix_obj{2};
        
    % calculate SVD coil combination
    prescanIm = prescan_twix_obj.image(:,:,:,:,1,1,1,1,1,1);
    prescanIm = ifft3s(permute(prescanIm,[4 3 1 2]));
    [~,~,V] =  svd(reshape(prescanIm,[],size(prescanIm,4)),'econ');
    reconPars.svdpars = V(:,1);
%     prescanIm_svd = reshape( reshape(prescanIm,[],size(prescanIm,4)) * reconPars.svdpars , size(prescanIm(:,:,:,1)) );
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
        timingReport{iAve} = reconstructSiemensVolumeTCL(twix_obj,thisReconPars);
    end
    
else
    if nRep > 1
        disp('Error - handling of data with multiple repetitions not yet implemented...! Please contact gallichand@cardiff.ac.uk for more info')
        return
    else
        timingReport = reconstructSiemensVolumeTCL(twix_obj,reconPars);       
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
fprintf('***** reconstructSiemensMP2RAGEwithTCL.m completed! *****\n')
fprintf('*************************************************************\n')
fprintf(['Total reconstruction time: ' num2str(totalTime_hrs) ' hours, ' num2str(totalTime_mins) ' mins\n']);

%%


    


end

