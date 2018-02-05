function [ACSims, timingReport, fatnavdir] = processFatNavs_GRAPPA4x4(twix_obj, outRoot, varargin)
% function [ACSims, timingReport, fatnavdir] = processFatNavs_GRAPPA4x4(twix_obj, outRoot, varargin)
%
% Process the FatNavs contained within the file specified by the 'twix_obj'
% (loaded in using 'mapVBVD_fatnavs.m') - using the folder 'outRoot' for
% output. A subfolder 'fatnavs_MIDXXX' will be created to store the
% individual FatNav images, and the resulting motion parameters saved as
% 'motion_parameters_spm_MIDXXX.mat' 
%
% If the output .mat file for the motion parameters already exists, it will
% skip all the processing and allow the existing file to be used.
%
% Options:-
%%%%%%%%%
%
%   'indexOfFirstNav' - In some data from testing various sequences the first
%                       FatNav was empty. This allows a later volume to be
%                       used as the reference for registration. Default is 1.
%
%   'showFigs'      - Choose whether or not to show progress to the screen
%                     in a figure window. Default is true.
%     
%   'nVirtualCoilsFatNavs' - decide how many coil compressed virtual coils
%                            to use for the FatNavs - default = all
%                            So far this doesn't really work for FatNavs :)
%  
%   'FatNavRes_mm' - able to specify other resolutions (but for now, must
%                    still be isotropic, 4x4 acceleration and same FOV -
%                    176x256x256mm - or 192x264x264mm for 6mm FatNavs)
%
%   'nSliceNeighbours' - include neighbouring virtual slices in determination 
%                        of GRAPPA weights. Some quick-and-dirty tests suggested 
%                        that '2' should be a good value to use for some standard 
%                        GRE sequences. For now, I'm leaving the default as
%                        '0' as it may affect registration results before
%                        and after this change...
%
%   'bApplyNoseFix' - Add option as to whether or not to attempt to
%                       down-weight the nose based on knowledge of
%                       locations of coil sensitivities for 32ch coil
%
%   'iAve'       - Which 'average' to reconstruct
%
%   'appendString' - String which gets appended to the output foldername
%
%  ------------
%  daniel.gallichan@epfl.ch, June 2015
%
%  Updates:
%
%  22/6/15 - danielg - changed name of variable 'indexOfRefNav' to
%         'indexOfFirstNav' and tried to fix handling of GIF creation when
%         this is being used.
%
%  4/1/16 - danielg - the NIFTI files for each FatNav now have the
%         resolution specified in them, which should make SPM realign behave more
%         reliably for lower resolutions. However, this has the effect of
%         CHANGING the definitions of motion parameters which are saved...!
%         Also added the new defaults for SPM realign of 'fwhm=3' and
%         'sep=2'.
% 
%  10/2/16 - danielg - added option to keep or delete FatNavs folder at the end
%              --> now moved this feature to part of
%              reconstructMP2RAGEiwthFatNavs.m instead - 15/3/16
%
%  8/3/16 - danielg - should now cope with 2, 4 or 6 mm FatNavs, so also
%         renamed the function to reflect this.
%
%  14/3/16 - danielg - now uses 'parfor' for each volume, and added
%          'nSliceNeighbours' option.
%
%  15/3/16 - danielg - make noseFix optional
%
%  21/4/16 - danielg - switch to parfor version of spm_realign
%                    - switch to SPM12
%                    - nosefix now default to off, as SPM12 allows wrap option
%
%  5/1/18 - gallichand@cardiff.ac.uk (danielg)
%         - returns the 'fatnavdir' folder
%         - incorporate 'average' counter to reconstruct only those FatNavs
%         -  NB: really the sequence itself needs updating on the scanner
%            so that the counters are incremented properly... Probably this
%            feature won't be used enough to make that worthwhile though...

if (~isfield(twix_obj,'FatNav') || ~isfield(twix_obj,'FatNav_refscan'))
    disp('Error: input twix_obj doesn''t seem to contain FatNav data...')
    return;
end


[indexOfFirstNav, bShowFigs, nVirtualCoilsFatNavs, FatNavRes_mm, nSliceNeighbours, bApplyNoseFix, iAve, appendString] = ...
    process_options(varargin,'indexOfFirstNav',1,'showFigs',1,'nVirtualCoilsFatNavs',[],...
                    'FatNavRes_mm',2,'nSliceNeighbours',0,'bApplyNoseFix',0, 'iAve', 1, 'appendString', '');

%% For testing code without calling as function:

% indexOfFirstNav = 1;
% bShowFigs = 1;
% nVirtualCoilsFatNavs = [];
% FatNavRes_mm = 4;
% nSliceNeighbours = 0;

%%

figIndex = 999; % use this figure number for live 'show' of how things are going...

Ax = 4; Ay = 4; % Hard-coded for a 'known' navigator...
nc_orig = twix_obj.image.NCha;
if ~isempty(nVirtualCoilsFatNavs)
    nc = nVirtualCoilsFatNavs;
else
    nc = nc_orig;
end
% FatNavRes_mm = 2; % For now this needs to be isotropic to keep things simple...

% make sure that oversampling is already removed here...
twix_obj.FatNav_refscan.flagRemoveOS = 1; 
twix_obj.FatNav.flagRemoveOS = 1;


% apply shift due to fat offset frequency (this assumes BW of 1950 Hz, and a fat/water shift of 3.5ppm)
offsetsXYZ = [0 0 (twix_obj.hdr.MeasYaps.sTXSPEC.asNucleusInfo{1,1}.lFrequency*3.5e-6/1950)]; 

nT = twix_obj.FatNav.dataSize(9)/twix_obj.FatNav.NAve;


%%

fileName = twix_obj.image.filename;
MIDstr = getMIDstr(fileName);

if outRoot(1)=='~'  % SPM seems to have a 'thing' about the tilde...
    outRoot = [getenv('HOME') outRoot(2:end)];
end

fatnavdir = [outRoot '/fatnavs_' MIDstr appendString];

if ~exist(fatnavdir,'dir')
    mkdir(fatnavdir);
end


 %%

dataSize = twix_obj.FatNav_refscan.dataSize;

switch FatNavRes_mm % in newer version of FatNav sequences, the resolution can be chosen at 2, 4 or 6 mm - with the FOV hard-coded in the sequence itself
    case {2,4}
        nx = 176/FatNavRes_mm;
        ny = 256/FatNavRes_mm;
        nz = 256/FatNavRes_mm;
    case 6
        nx = 192/FatNavRes_mm;
        ny = 264/FatNavRes_mm;
%         nz = 264/FatNavRes_mm;
        nz = 64; % No idea why, but the 6mm data has 64 points in the readout direction for the ACS lines instead of 44...
end

ACSdata = twix_obj.FatNav_refscan();
ACSdata = squeeze(ACSdata);
ACSdata = ACSdata(:,:,:,:,end); % from some old data the FatNavs might have a 'set' counter of 2
if nc_orig ~= nc       
    % calculate coil compression matrices
    % using code from Miki Lustig's GCC package
    % (http://www.eecs.berkeley.edu/~mlustig/Software.html)
    % which implements the method from Zhang et. al MRM 2013;69(2):571-82.
    ACSdata = permute(ifft1s(ACSdata,1),[1 3 4 2]);
    nread = size(ACSdata,1); npe1 = size(ACSdata,2); npe2 = size(ACSdata,3);
    
    mtx = zeros(nc_orig,min(nc_orig,npe1*npe2),nread);

    for iRead = 1:nread
        iiRead = iRead-4:iRead+4; % use neighbouring points in readout direction to make weights smoother
        iiRead(iiRead<1) = [];
        iiRead(iiRead>nread) = [];
        tmpc = reshape(ACSdata(iiRead,:,:,:),length(iiRead)*npe1*npe2,nc_orig);
        [~,~,V] = svd(tmpc,'econ');
        mtx(:,:,iRead) = V;
    end
    
     mtx = mtx(:,1:nVirtualCoilsFatNavs,:); % cut down to just number of virtual coils
    
    % now align those compression matrices (from alignCCMtx.m)
    
    % align everything based on the middle slice.
    n0 = floor(nread/2);
    A00 = mtx(:,1:nVirtualCoilsFatNavs,n0);
    
    % Align backwards to first slice
    A0 = A00;
    for iRead = n0-1:-1:1
        A1 = mtx(:,1:nVirtualCoilsFatNavs,nread);
        C = A1'*A0;
        [U,~,V]= svd(C,'econ');
        P = V*U';
        mtx(:,1:nVirtualCoilsFatNavs,nread) = A1*P';
        A0 = mtx(:,1:nVirtualCoilsFatNavs,nread);
    end
    
    % Align forward to end slice
    A0 = A00;
    for n = n0+1:nread
        A1 = mtx(:,1:nVirtualCoilsFatNavs,nread);
        C = A1'*A0;
        [U,~,V]= svd(C,'econ');
        P = V*U';
        mtx(:,1:nVirtualCoilsFatNavs,nread) = A1*P';
        A0 = mtx(:,1:nVirtualCoilsFatNavs,nread);
    end
   
    % and directly apply to the ACS data:
    newACSdata = zeros(nread,npe1,npe2,nVirtualCoilsFatNavs);
    for iRead = 1:nread
        thistemp = reshape(ACSdata(iRead,:,:,:),npe1*npe2,nc_orig)*mtx(:,:,iRead);
        newACSdata(iRead,:,:,:) = reshape(thistemp,[1 npe1 npe2 nVirtualCoilsFatNavs]);
    end
    ACSdata = fft1s(reshape(permute(newACSdata,[1 4 2 3]),[nread nVirtualCoilsFatNavs npe1 npe2]),1);       

    % now interpolate weights as for some reason ACS data for FatNavs has
    % only half the number of readout points...
    mtx = reshape(mtx,nc_orig*nVirtualCoilsFatNavs,nz/2);
    newmtx = zeros(nc_orig*nVirtualCoilsFatNavs,nz);
    for iM = 1:nc_orig*nVirtualCoilsFatNavs
        newmtx(iM,:) = interp1(linspace(0,1,nz/2),mtx(iM,:),linspace(0,1,nz));
    end
    mtx = reshape(newmtx,nc_orig,nVirtualCoilsFatNavs,nz);
else
    mtx = []; % parfor needs this to exist to be happy...
end
    
    
    
ACSdata = permute(ACSdata,[4 3 1 2]);
ACSdata(nx,ny,nz,1) = 0; % extend to full size
ACSdata = circshift(ACSdata,double(round([nx/2-dataSize(4)/2 ny/2-dataSize(3)/2 0 0])));  
% no idea why it should be necessary to put 'double' here - but I got the error 'invalid shift type: must be a finite, nonsparse, real integer vector

[ox, oy, oz] = ndgrid(linspace(-1,1,nx),linspace(-1,1,ny),linspace(-1,1,nz));
offsetPhase = exp(1i*pi*( offsetsXYZ(1)*ox + offsetsXYZ(2)*oy + offsetsXYZ(3)*oz) );
ACSdata = bsxfun(@times,ACSdata,offsetPhase);

ACSdata = ifft1s(ACSdata,3);
ACSdata = ACSdata(:,:,end:-1:1,:);

ACSims = ifft2s(ACSdata);
ACSim = ssos(ACSims);

ovACS = orthoview(ACSim,'drawIms',0);
imab_overwrite([fatnavdir '/a_FatNav_ACSim_' MIDstr '.png'],ovACS.oneIm);

ACSdata = zeroCrop(ACSdata); % after making an image at full resolution, drop back to just where the data really is

if FatNavRes_mm==2 && bApplyNoseFix
    % only attempt attenuation of aliased nose signal for 2mm data...
    if strcmp(twix_obj.hdr.MeasYaps.asCoilSelectMeas{1}.asList{1}.sCoilElementID.tCoilID,'32Ch_Head_7T')
        disp('...............')
        disp('Detected use of 32Ch_Head_7T RF coil: ')
        disp('Will also attempt to use prior knowledge of the locations of the channels in order to down-weight')
        disp('signal from the nose')
        disp('...............')
        %     iC_keep = [     4     6    10    14    20    22    27    31     1     5    13    21    28    32     2     9    17    18     3     7    11    15    19    23    26    30];
        iC_keep = 1:nc; % was previously using 26 upper channels, but recon will be faster to use fewer compressed virtual coils
        % If really desired, could think about putting in
        % knowledge of coil positions into chosen virtual coils
        
        %%% sort out the nose...
        coilMaxCoords =    [31    98    83;    27    86    87;    27   111    62;    37   113    68;    24    86    85;     8    75    62;...
            14    89    48;     9    71    26;    40    25    80;    37    24    76;    24    23    54;    29    22    34;    20    38    80;...
            14    38    67;    10    43    49;     7    64    27;    44    83    94;    57   103    81;    62   112    61;    53   113    69;...
            68    82    86;    78    90    65;    78    88    49;    80    95    30;    84    62    25;    83    51    54;    82    56    67;...
            60    58    92;    59    21    34;    66    23    53;    62    25    71;    64    33    82].';
        
        [xx, yy, zz] = ndgrid(1:nx,1:ny,1:nz);
        sigDist = zeros(size(ACSims));
        for iC = 1:nc
            this_r = sqrt( (xx-coilMaxCoords(1,iC)).^2 + (yy-coilMaxCoords(2,iC)).^2 + (zz-coilMaxCoords(3,iC)).^2 );
            sigDist(:,:,:,iC) = (this_r.^3).*abs(ACSims(:,:,:,iC));
        end
        sigDistMask = zeros(size(sigDist));
        sigDistMask(sigDist<(150*norm(ACSims(:)))) = 1;
        
        ov1 = orthoview(ssos(1-sigDistMask),'mip',1,'drawIms',0);
        imab_overwrite([fatnavdir '/a_noseRemovalMask_mip.png'],ov1.oneIm);
        
    else
        iC_keep = 1:nc;
        disp('...............')
        disp(['RF coil unrecognized, using data from all ' num2str(nc) ' channels.'])
        disp(['No information available to reduce signal from nose'])
        disp('...............')
        sigDistMask = [];
        
    end
    
else
    iC_keep = 1:nc;
    sigDistMask = [];
end

%% Need to account for the fact that the resolution doesn't match for the 
%% 6mm data in the readout direction between the ACS data and the actual data
%% -- try to fix with interpolation in image space...

if FatNavRes_mm==6
    nz = 44;
    ACSims = ifft2s(ACSdata);
    szACS = size(ACSims);
    ACSmat = reshape(permute(ACSims,[3 1 2 4]),64,[]);
    ACSmat = interp1(linspace(0,1,64),ACSmat,linspace(0,1,nz));
    ACSims = permute(reshape(ACSmat,[nz szACS([1 2 4])]),[2 3 1 4]);
    ACSdata = fft2s(ACSims);
    
    offsetPhase = reshape(interp1(linspace(0,1,64),reshape(offsetPhase,[],64).',linspace(0,1,nz)).',nx,ny,nz);
end
    
    
%% Check if FatNavs have already been processed, and use them instead if they have

if exist([outRoot '/motion_parameters_spm_' MIDstr appendString '.mat'],'file');
    disp(['Existing estimated motion-parameters found - using them instead of recalculating...'])
    timingReport.calculateGRAPPAweights = 0;
    timingReport.eachFatNav = zeros(nT,1);
    timingReport.SPMalignment = 0;
    timingReport.allFatNavs = 0;
    return
end

nProcessedImages = length(dir([fatnavdir '/eachFatNav_*.nii']));

if nProcessedImages>0
    disp(['Found ' num2str(nProcessedImages) ' existing processed FatNavs and expected to find ' num2str(nT)])
    timingReport.calculateGRAPPAweights = 0;
    timingReport.eachFatNav = zeros(nT,1);
    timingReport.SPMalignment = 0;
    timingReport.allFatNavs = 0;
end

if nProcessedImages~=nT % assume images are there, but not the alignment of them
    
    %% calculate and save the GRAPPA weights
    
    gx = 2; gy = 2; % size of GRAPPA kernel
    grapKernel = zeros(gx + (gx-1)*(Ax-1), gy + (gy-1) * (Ay-1));
    
    % source points are marked with 1
    grapKernel(1:Ax:end, 1:Ay:end) = 1;
    
    % target points are marked with 0.5
    if gx == 1,    startx = 2; else  startx = 2 + (floor(gx/2) - 1) * Ax; end
    if gy == 1,    starty = 2; else  starty = 2 + (floor(gy/2) - 1) * Ay; end
    
    grapKernel(startx:startx+Ax-2, starty:starty+Ay-2) = 0.5;
    grapKernel(startx:startx+Ax-1, starty:starty+Ay-2) = 0.5;
    grapKernel(startx:startx+Ax-2, starty:starty+Ay-1) = 0.5;
    
    stepSize = [Ax Ay];
    all_grapW = zeros(gx*gy*nc,((Ax*Ay)-1)*nc,128);
    zMax = nz; % don't need to go outside the brain... would save time to detect the top of the brain here - but could be difficult to make robust...
    
    tic
    disp('... Calculating GRAPPA weights ....')
    parfor iZ = 1:zMax
        
        if nSliceNeighbours > 0
            
            src = []; targ = [];
            
            for this_iZ = iZ-nSliceNeighbours:iZ+nSliceNeighbours
                if this_iZ > 1 && this_iZ <= nz
                    thisACS = squeeze(ACSdata(:,:,this_iZ,:));
                    [~, this_src, this_targ] = GrappaCalib3D_arb(thisACS,grapKernel,0,0);
                    src = [src; this_src];
                    targ = [targ; this_targ];
                end
            end
            
            all_grapW(:,:,iZ) = pinv(src)*targ;
            
        else
            thisACS = squeeze(ACSdata(:,:,iZ,:));
            all_grapW(:,:,iZ) = GrappaCalib3D_arb(thisACS,grapKernel,0,1);
        end
        
        fprintf('.');
    end
    fprintf('\n');
    disp('...Done...')
    timingReport.calculateGRAPPAweights = toc;
    
    
    %%
    
    nxMeas = twix_obj.FatNav.dataSize(4);
    nyMeas = twix_obj.FatNav.dataSize(3);
    startPos = [nx-nxMeas+1 ny-nyMeas+1];
    
    timingEachFatNav  = zeros(nT,1);
    timeBeforeFatNavs = clock;
    
    parfor iT = indexOfFirstNav:nT
        
        full_GRAPPArecons = int16(zeros(nx,ny,nz));
        
        fatnavdata = twix_obj.FatNav(:,:,:,:,1,iAve,1,1,iT+(iAve-1)*nT);
        if nc~=nc_orig
            % do coil compression
            fatnavdata = permute(ifft1s(fatnavdata,1),[1 3 4 2]);
            [nread,npe1,npe2] = size(fatnavdata(:,:,:,1));
            newfatnavdata = zeros(nread,npe1,npe2,nVirtualCoilsFatNavs);
            for iRead = 1:nread
                thistemp = reshape(fatnavdata(iRead,:,:,:),npe1*npe2,nc_orig)*mtx(:,:,iRead);
                newfatnavdata(iRead,:,:,:) = reshape(thistemp,[1 npe1 npe2 nVirtualCoilsFatNavs]);
            end
            fatnavdata = fft1s(reshape(permute(newfatnavdata,[1 4 2 3]),[nread nVirtualCoilsFatNavs npe1 npe2]),1);
        end
        
        fatnavdata = permute(fatnavdata,[4 3 1 2]);
        fatnavdata(nx,ny,nz,1) = 0; % extend to full size
        
        fatnavdata = circshift(fatnavdata,double([nx-nxMeas ny-nyMeas 0 0]));
        
        if any(offsetsXYZ)
            fatnavdata = bsxfun(@times,fatnavdata,offsetPhase);
        end
        
        fatnavdata = ifft1s(fatnavdata,3);
        fatnavdata = fatnavdata(:,:,end:-1:1,:);
        
        tic
        for iZ = 1:zMax
            thisfatnavdata = squeeze(fatnavdata(:,:,iZ,:));
            thisfatnavdata(nx,ny,1) = 0;
            
            grapK = GrappaReco3D_arb(thisfatnavdata,grapKernel,all_grapW(:,:,iZ),stepSize,startPos);
            
            f1 = tukey(nx,.7,.3); f2 = tukey(ny,.8,.3);
            grapK = bsxfun(@times,grapK,f1.');
            grapK = bsxfun(@times,grapK,f2);
            
            grapIms = ifft2s(grapK);
            if isempty(sigDistMask)
                grapIm = ssos(grapIms(:,:,iC_keep))*nx*ny*nz;
            else
                grapIm = ssos(grapIms(:,:,iC_keep).*squeeze(sigDistMask(:,:,iZ,iC_keep)))*nx*ny*nz;
            end
            full_GRAPPArecons(:,:,iZ) = int16(4095*grapIm/20); % use 20 as first signal to clip...
            
        end
        disp(['FatNav ' num2str(iT) ' out of ' num2str(nT)]);
        timingEachFatNav(iT) = toc;
        
        
        sn(full_GRAPPArecons,[fatnavdir '/eachFatNav_' num2str(iT,'%.3d') '.nii'],FatNavRes_mm*[1 1 1]);
        
        if bShowFigs
            fig(figIndex)
            orthoview(full_GRAPPArecons,'useNewFig',0);
            set(gcf,'Position',[    50   720   950  340]);
            subplot1(1)
            ylabel(['FatNav ' num2str(iT) ' out of ' num2str(nT)]);
            drawnow;
        end
        
    end
    
    timeAfterFatNavs = clock;
    timingReport.allFatNavs = etime(timeAfterFatNavs,timeBeforeFatNavs);
    
    timingReport.eachFatNav = timingEachFatNav;
    
end

%% Do registration
   

% do registration with SPM

if indexOfFirstNav>1
    iiV = indexOfFirstNav:nT;
    nT = nT-indexOfFirstNav+1;
else
    iiV = 1:nT;
end

for iV = 1:nT
    V(iV) = spm_vol_nifti([fatnavdir '/eachFatNav_' num2str(iiV(iV),'%.3d') '.nii']);
end


alignpars = struct('quality',1,'fwhm',3,'sep',2,'rtm',0,'wrap',[0 1 0]);
disp('...............')
disp(['Aligning ' num2str(nT) ' FatNavs using SPM'])
tic
% PPos = spm_realign(V,alignpars);
PPos = spm12_realign_parfor(V,alignpars);
disp(['Done!'])
timingReport.SPMalignment = toc;
disp('...............')
MPos = get_spm_affmats(PPos);

fitpars = MPos.pars(1:6,:);
% fitpars(1:3,:) = fitpars(1:3,:)*FatNavRes_mm; % put into mm <-- now already in mm
fitpars(4:6,:) = fitpars(4:6,:)*180/pi;

MPos_cent = get_spm_affmats(PPos,round(nT/2+1));

fitpars_cent = MPos_cent.pars(1:6,:);
% fitpars_cent(1:3,:) = fitpars_cent(1:3,:)*FatNavRes_mm; % put into mm <-- now already in mm
fitpars_cent(4:6,:) = fitpars_cent(4:6,:)*180/pi;

save([outRoot '/motion_parameters_spm_' MIDstr appendString '.mat'],'MPos','fitpars','MPos_cent','fitpars_cent');
            
     
%%

fig(figIndex);
plotFitPars(fitpars_cent);
export_fig([fatnavdir '/a_FatNav_MoCoPars_' MIDstr '.png'])


%% make gif movies of before and after alignment (requires ImageMagick to make the GIF)

disp('...............')
disp('... Checking for ImageMagick...')
testMagick = system('convert -version');
if testMagick==0 % apparently a status of '0' means that everything is good...
    disp('...............')
    disp('... ImageMagick found, you will get some animated gifs too...')
    disp('...............')
    iOut = unique(round(linspace(1,nT,15))); % up to 15 images, if we have enough
    
    VfromPos = apply_spm_affmats(V,MPos.mats);
    spm_reslice(VfromPos(iOut),struct('mask',false,'mean',false,'interp',5,'which',2,'wrap',[0 0 0],'prefix','spm_'));
    
    for iiOut = 1:length(iOut)
        thisOutIm = rn([fatnavdir '/eachFatNav_' num2str(iOut(iiOut)+(indexOfFirstNav-1),'%.3d') '.nii']);
        if iiOut==1
            climsMax = percentile(thisOutIm(:),99);
        end
        ov1 = orthoview(thisOutIm,'drawIms',0);
        imab_overwrite([fatnavdir '/temp_mov_' num2str(iiOut,'%.3d') '.png'],ov1.oneIm,[0 climsMax]);
    end   
    processString = ['convert -dispose 2 -delay 10 -loop 0 ' fatnavdir '/temp_mov_*.png ' fatnavdir '/a_mov_eachFatNav.gif'];
    system(processString);
    delete([fatnavdir '/temp_mov_*.png'])
    
    for iiOut = 1:length(iOut)        
        thisOutIm = rn([fatnavdir '/spm_eachFatNav_' num2str(iOut(iiOut)+(indexOfFirstNav-1),'%.3d') '.nii']);
        ov1 = orthoview(thisOutIm,'drawIms',0);
        imab_overwrite([fatnavdir '/temp_mov_' num2str(iiOut,'%.3d') '.png'],ov1.oneIm,[0 climsMax]);
    end
    processString = ['convert -dispose 2 -delay 10 -loop 0 ' fatnavdir '/temp_mov_*.png ' fatnavdir '/a_mov_spm_eachFatNav.gif'];
    system(processString);    
    delete([fatnavdir '/temp_mov_*.png'])
    
else
    disp('...............')
    disp('... ImageMagick not found, install this from www.imagemagick.org and re-run if you would like animated gifs')
    disp('...............')
end

