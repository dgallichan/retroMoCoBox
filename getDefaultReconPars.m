function reconPars = getDefaultReconPars(computeRAMlevel,B0_fieldStrength_Tesla, sequenceType)
% function reconPars = getDefaultReconPars(computeRAMlevel,B0_fieldStrength_Tesla, sequenceType)
% 
% Gets the default set of reconstruction parameters

reconPars.rawDataFile = ''; 
% 'rawDataFile' - this is the filename for the raw data you exported from
%                 TWIX (including the full path). Due to the way I chose to
%                 name output files, if you have renamed the datafile it is
%                 required to still contain the text 'MID' followed by a
%                 number and then an underscore ('_') so that multiple
%                 files can be processed from the same directory and not
%                 get overwritten.

reconPars.retroMocoBoxVersion = retroMocoBoxVersion(); % for tracking version

if nargin < 1 || strcmp(computeRAMlevel,'normalRAM')
    
    reconPars.bGRAPPAinRAM = 0;
%      'bGRAPPAinRAM' - set this to '1' to perform the GRAPPA recon (if
%                       necessary!) for the host sequence entirely in RAM. This
%                       can be quite a lot faster - but requires sufficient
%                       RAM...! 
%                       (Note that the temporary variables for applying the
%                       MoCo are currently done outside of RAM as this
%                       is currently not the bottleneck for that part of
%                       the recon)

    reconPars.bKeepReconInRAM = 0;
%      'bKeepReconInRAM' - set this to '1' to keep all the reconstructed
%                          files (INV1, INV2 and UNI, all with and without
%                          correction) in RAM, or default is to use a
%                          temporary file.

    reconPars.bFullParforRecon = 0;
%      'bFullParforRecon' - set this to '1' to enable parfor over the NUFFT
%                           loop. For this to work, bKeepReconInRAM must
%                           also be set - and depending how many CPUs you
%                           have available on your matlabpool, you may need
%                           to have LOADS of RAM...  But it is much faster!

elseif strcmp(computeRAMlevel, 'highRAM')
    reconPars.bGRAPPAinRAM = 1;
    reconPars.bKeepReconInRAM = 1;
    reconPars.bFullParforRecon = 0;

elseif strcmp(computeRAMlevel, 'veryhighRAM')
    reconPars.bGRAPPAinRAM = 1;
    reconPars.bKeepReconInRAM = 1;
    reconPars.bFullParforRecon = 1;
else
    error('Invalid computeRAMlevel choice')
end
   


if nargin < 2
    B0_fieldStrength_Tesla = []; % don't assume a field strength
    % (but, this will still currently assume a FatNav resolution of 2mm later on...)
    %
    % Used for assumptions:
    % 7T -> FatNav Res 2mm for MPRAGE, 4mm for SPACE or GRE
    % 3T -> FatNav Res 4mm for all
end

if nargin < 3
    sequenceType = 'MPRAGE';
end





reconPars.outRoot = [];
% 'outRoot' -    this gets used as the output folder for all processing
%                 (and the default root for the temporary folder which may
%                 generate several Gb of temporary files along the way - 
%                 - override this with 'tempRoot' option). The default 
%                 'outRoot' is to use the same folder as the raw data file,
%                 achieved by using '[]';

reconPars.outFolderPrefix = sequenceType;
%      'outFolderPrefix' - what to put in front of the MIDstr in the
%                          generated output folder - default is the
%                          sequence type

reconPars.tempRoot = [];
%      'tempRoot' - location for creating temporary files - can be many Gb.
%                   Default is to use the same location as 'outRoot',
%                   achieved by using '[]';


reconPars.bLinParSwap = 0;
%      'bLinParSwap' - set this to '1' to indicate that the 'LIN/PAR swap'
%                     option was chosen in the MP2RAGE sequence. This
%                     alters which direction in k-space the FatNavs
%                     correspond to.


reconPars.bKeepGRAPPArecon = 0;
%      'bKeepGRAPPArecon' - set this to '1' to prevent deleting of the
%                           GRAPPA recon of the host data (which would be
%                           deleted by default). 

reconPars.coilCombineMethod = 'default';
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


if B0_fieldStrength_Tesla < 4
    reconPars.FatNavRes_mm = 4;
elseif strcmp(sequenceType,'MPRAGE')
    reconPars.FatNavRes_mm = 2;
else
    reconPars.FatNavRes_mm = [];
end
%       'FatNavRes_mm' - the spatial resolution of the acquired FatNavs, 
%                        specified in mm. The current implementation of the
%                        pulse sequence does not allow this information to
%                        be stored in the raw data, so this must be entered
%                        here manually - or use the default of 2 mm.
%                        Whatever you specify here, the code does still
%                        assume that the acceleration for the FatNavs was
%                        4x4 and that the FOV was 176x256x256mm for 2mm and
%                        4mm (and slightly different for 6mm)

reconPars.flipAxes_xyz = [];
%       'flipAxes_xyz' - specify which axes should be flipped for the final
%                        NIFTI outputs. Note that unlike the previous way
%                        that 'swapDims_xyz' worked, the flipAxes_xyz has
%                        no effect on aligning FatNavs with the host
%                        sequence, and is purely for 'aesthetics' of how
%                        the final 3D volumes look (MoCo is now applied in
%                        the 'native' RPS coordinates of the host sequence,
%                        before transforming to XYZ). Default is [] which
%                        uses default of 'reconstructSiemensVolume.m' which
%                        is currently [0 0 1].

reconPars.bZipNIFTIs = 1;
%       'bZipNIFTIs' - Use '1' to apply gzip at the end to all the NIFTI
%                     files (default), otherwise just leave them uncompressed.

reconPars.bKeepFatNavs = 0;
%       'bKeepFatNavs' - Use '1' to keep all the reconstructed FatNavs,
%                       otherwise delete that folder when finished
%                       (default).

reconPars.bKeepComplexImageData = 0;
%       'bKeepComplexImageData' - save out the complex data per coil as MATLAB files 
%                                 before and after application of MoCo.

reconPars.bKeepPatientInfo = 1;
%       'bKeepPatientInfo' - Use '1' to keep sensitive info from raw data header  
%                           in HTML output (default). Use the number 0 to anomyize 
%                           completely or use a string e.g. '0019' to insert the ID from
%                           another database.

reconPars.bUseGPU = 0;
%       'bUseGPU' - decides whether to use gpuNUFFT
%                   (from https://github.com/andyschwarzl/gpuNUFFT and needs 
%                   to be on MATLAB path) or default CPU using MIRT NUFFT

reconPars.parpoolSize = feature('numcores'); 
%        'parpoolSize' - on multi-CPU codes can explicitly set the parpool
%                        size if desired (or start manually in advance).

reconPars.NUFFTosf = 1.5;
%        'NUFFTosf' - the oversampling factor in the NUFFT. Larger numbers
%        are slower and more RAM intensive, but potentially more accurate.
%        A value of 1.5 seems a good compromise in my experience

reconPars.bSwapFatNavHandedness = 1;
%        'bSwapFatNavHandedness' - because the FatNav is acquired without
%        any slice prescription from the scanner operator, it seems that it
%        may still need a L/R flip for natural alignment against the host
%        sequence. This is currently enabled by default.

reconPars.bApplyFatNavNoseCircshift = 1;
%   'bApplyFatNavNoseCircshift' - Finds the minimum signal in FatNav in PE
%                           direction and circularly shifts by this amount
%                           to attempt to centre the FOV for more reliable
%                           motion estimates. It then (hopefully!) accounts
%                           for this offset properly in the motion
%                           estimates, keeping them relative to isocentre.

reconPars.bApplyFatNavNeckMasking = 1;
%   'bApplyFatNavNeckMasking' - Aligns a FatNav to a 'standard' FatNav (included in
%                       'fatnavStandardSpace' folder) and puts a mask
%                       defined in that space back into subject space to
%                       try to get alignment to concentrate only on region
%                       of image where FatNav can be expected to be assumed
%                       rigid.

reconPars.bGetGRAPPA_SVD = 0;
%    'bGetGRAPPA_SVD' - Useful mainly for debugging for GRAPPA to also return a 
%                       single virtual coil image with signal almost everywhere
    
reconPars.extraFlipMat1 = diag([-1 1 -1]);
    % This is an additional flip of axes compared to what comes out of the
    % ported Siemens code. Seems there are always reasons this might
    % genuinely be necessary, so hopefully this can stay as diag([-1 1 -1)
    % for all data.

reconPars.extraFlipMat2 = eye(4);
    % This is a 'just in case' matrix that allows further flips or
    % manipulations of the motion parameters themselves. Hopefully can stay
    % as eye(4) for all data.

reconPars.extraPositionOffsetSignFlips = [1 -1 1];
    % This is because I'm *still* having trouble getting FatNavs to align
    % with host sequence if the host sequence is rotated e.g. PA instead of AP

