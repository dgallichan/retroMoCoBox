function timingReport = reconstructSiemensMPRAGEvolume(twix_obj,reconPars)
% function timingReport = reconstructSiemensMPRAGEvolume(twix_obj,reconPars)
% 
% Called by reconstructSiemensMP2RAGEwithFatNavs.m to reconstruct a
% specific 'average' or 'repetition'.
%
%
% -- July 2026, gallichand@cardiff.ac.uk -> major overhaul for v1.0.0dev
%     see retroMocoBoxVersion.m for more info

retroMocoBoxVersion = reconPars.retroMocoBoxVersion; % put this into the HTML for reference

%%

if ~isfield(reconPars,'iAve')
    reconPars.iAve = 1;
end
if ~isfield(reconPars,'iRep')
    reconPars.iRep = 1;
end
if isfield(reconPars,'swapDims_xyz')
    disp("WARNING: swapDims_xyz is deprecated, use flipAxes_xyz instead - name has change as functionality has changed")
end

if ~isfield(reconPars,'flipAxes_xyz') || isempty(reconPars.flipAxes_xyz)
    % Previously the 'swapDims_xyz' also affected the MoCo, the new
    % 'flipAxes_xyz' is purely aesthetic for what space the output volume
    % is put into.
    % flipAxes_xyz = [1 1 0]; % I currently expect this to be correct, as z is oriented 'upside down' for Siemens, and they use FFT instead of iFFT
    reconPars.flipAxes_xyz = [0 0 1]; % this is looking better so far!
end

outputSuffix = 'MoCo';

%%           

%%% Make local function as inline so that it can be used in scripts
local_reformatDateString = @(dateString) [dateString(7:8) '/' dateString(5:6) '/' dateString(1:4)];

%%

[datadir, fileName, fileExt] = fileparts(reconPars.rawDataFile);
if isempty(reconPars.outRoot)
    reconPars.outRoot = datadir;
end
if reconPars.outRoot(1)=='~'  % SPM seems to have a 'thing' about the tilde...
    reconPars.outRoot = [getenv('HOME') reconPars.outRoot(2:end)];
end

MIDstr = getMIDstr(reconPars.rawDataFile);

filterFrac = 0.05; % used when 'lowres' coil combination method is selected (not default..)

startTime = clock;

%% Make output folders

appendString = '';
if reconPars.iAve > 1
    appendString = ['_Ave' num2str(reconPars.iAve)];
end
if reconPars.iRep > 1
    appendString = [appendString '_Rep' num2str(reconPars.iRep)];
end

outDir = [reconPars.outRoot '/' reconPars.outFolderPrefix '_' MIDstr appendString];
if ~exist(outDir,'dir')
    mkdir(outDir)
end

htmlDir = [outDir '/html'];
if ~exist(htmlDir,'dir')
    mkdir(htmlDir)
end

if isempty(reconPars.tempRoot)
    reconPars.tempRoot = outDir;
end
tempDir = [reconPars.tempRoot '/temp_' MIDstr appendString];
if ~exist(tempDir,'dir')
    mkdir(tempDir)
end

% intialize html index file
fidHTML = fopen([htmlDir '/index.html'],'w');
fprintf(fidHTML,['<html><head><title>MP(2)RAGE with FatNavs - Summary</title>\n']);
fprintf(fidHTML,'</head>\n');
fprintf(fidHTML,'<body>\n');
fprintf(fidHTML,['<h2>MP(2)RAGE with FatNavs - Summary: %s</h2>\n'],[fileName]);

% include version number
fprintf(fidHTML,['<em>' char(datetime) '- created with reconstructSiemensMP2RAGEwithFatNavs.m, version: ' retroMocoBoxVersion '</em><br><br>\n']);






%%

if isempty(reconPars.FatNavRes_mm)
    manualFatNavRes = 0;
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
    manualFatNavRes = 1;
end

switch reconPars.FatNavRes_mm % in newer version of FatNav sequences, the resolution can be chosen at 2, 4 or 6 mm - with the FOV hard-coded in the sequence itself
    % these choices should also match the corresponding code in processFatNavs_GRAPPA4x4.m
    case {2,4}
        FatNav_FOVxyz = [176 256 256]; % FatNav FOV         
        FatNavDims_xyz = FatNav_FOVxyz ./ reconPars.FatNavRes_mm;
    case 6
        FatNav_FOVxyz = [192 264 264]; % FatNav FOV         
        FatNavDims_xyz = FatNav_FOVxyz ./ reconPars.FatNavRes_mm;
        FatNavDims_xyz(3) = 64; % No idea why, but the 6mm data has 64 points in the readout direction for the ACS lines instead of 44...        
end



%%


% in the next line I assume that Siemens always use the same structure for
% the 'tReferenceImage0' field - but I haven't looked for documentation to
% support this, so it may not always extract the scan date properly...
if isfield(twix_obj.hdr.MeasYaps,'tReferenceImage0') % some files apparently might not even have this field...
    iDot = strfind(char(twix_obj.hdr.MeasYaps.tReferenceImage0),'.'); % some versions seem to create a char already, others not...
    thisDateString = char(twix_obj.hdr.MeasYaps.tReferenceImage0);
    try
        if strcmp(thisDateString(iDot(end)+(1:6)),'300000') % it seems sometimes this 300000 is present, then the year is only two digits...
            fprintf(fidHTML,['<strong>Date of scan:</strong> ' local_reformatDateString(['20' thisDateString(iDot(end)+6+(1:6))]) '<br>\n']);
        else
            fprintf(fidHTML,['<strong>Date of scan:</strong> ' local_reformatDateString(thisDateString(iDot(end)+(1:8))) '<br>\n']);
        end
    catch
        fprintf(fidHTML,['<strong>Date of scan:</strong> unable to read from file <br>\n']);
    end
end

if ischar(reconPars.bKeepPatientInfo)
    fprintf(fidHTML,['<strong>Database ID:</strong> ' reconPars.bKeepPatientInfo '<br>\n']);
else
    if reconPars.bKeepPatientInfo
        %% Put patient info into HTML
        
      if isfield(twix_obj.hdr.Config,'PatientName')
            fprintf(fidHTML,['<strong>Patient Name:</strong> ' twix_obj.hdr.Config.PatientName '<br>\n']);
        else
            if isfield(twix_obj.hdr.Config,'tPatientName') % why would this be different for some scanners / software versions...?!?
                fprintf(fidHTML,['<strong>Patient Name:</strong> ' twix_obj.hdr.Config.tPatientName '<br>\n']);
            end
        end
        fprintf(fidHTML,['<strong>Patient ID:</strong> ' twix_obj.hdr.Config.PatientID '<br>\n']);
        switch twix_obj.hdr.Config.PatientSex
            case 1
                fprintf(fidHTML,['<strong>Patient Sex:</strong> Female<br>\n']);
            case 2
                fprintf(fidHTML,['<strong>Patient Sex:</strong> Male<br>\n']);
            case 3
                fprintf(fidHTML,['<strong>Patient Sex:</strong> Other<br>\n']);
        end
        try
            fprintf(fidHTML,['<strong>Patient date of birth:</strong> ' local_reformatDateString(num2str(twix_obj.hdr.Config.PatientBirthDay)) '<br>\n']);
        catch
            fprintf(fidHTML,['<strong>Patient date of birth:</strong> ' twix_obj.hdr.Config.PatientBirthDay '<br>\n']);
        end
        
    else
        fprintf(fidHTML,['<strong>Patient Info:</strong> Anonymised<br>\n']);
    end
end


%% Try to establish slice orientation - this has not been rigourously tested for all possible orientations...!


rotAndShift = getSiemensRotMatAndShift(twix_obj.hdr);

switch rotAndShift.baseOrientation 
    case 'Coronal'
        % danielg -  11/6/26 - added 'abs' and  >135 here to try to handle 180 deg in-plane rot!
        if abs(rotAndShift.dInPlaneRot*180/pi) < 45 || abs(rotAndShift.dInPlaneRot*180/pi) > 135
            permutedims = [2 3 1];
        else
            permutedims = [1 3 2]; 
        end       
    case 'Transverse'
        if abs(rotAndShift.dInPlaneRot*180/pi) < 45 || abs(rotAndShift.dInPlaneRot*180/pi) > 135
            permutedims = [1 2 3];
        else
            permutedims = [2 1 3];
        end
    case 'Sagittal'
        if abs(rotAndShift.dInPlaneRot*180/pi) < 45 || abs(rotAndShift.dInPlaneRot*180/pi) > 135
            permutedims = [3 2 1];
        else
            permutedims = [3 1 2];
        end
end


[~, permutedims_toXYZ] = sort(permutedims);

% Note that in the following Hrps refers to the dimensions of the host
% sequence when fully reconstructed, whereas hrps refers to the dimensions
% of the acquired data. To avoid wasting memory, matrices are not always
% filled to their full size immediately.
Hrps = [twix_obj.hdr.MeasYaps.sKSpace.lBaseResolution twix_obj.hdr.MeasYaps.sKSpace.lPhaseEncodingLines twix_obj.hdr.MeasYaps.sKSpace.lPartitions].';
% the line above can be adapted to get the right numbers if there is also oversampling in the slice direction:
% Hrps = [twix_obj.hdr.MeasYaps.sKSpace.lBaseResolution twix_obj.hdr.MeasYaps.sKSpace.lPhaseEncodingLines round(twix_obj.hdr.MeasYaps.sKSpace.lPartitions/(1+twix_obj.hdr.MeasYaps.sKSpace.dSliceOversamplingForDialog))].';
% but then other code would need to be changed to handle this properly...
% For the time-being I will recommend people not to use slice oversampling!
Hxyz = Hrps(permutedims); % full reconstructed matrix size

hrps = [twix_obj.image.NCol/2 twix_obj.image.NLin twix_obj.image.NPar].';
hxyz = hrps(permutedims);

Arps = [1 twix_obj.hdr.MeasYaps.sPat.lAccelFactPE twix_obj.hdr.MeasYaps.sPat.lAccelFact3D].';
Axyz = Arps(permutedims);

FOVrps = [twix_obj.hdr.MeasYaps.sSliceArray.asSlice{1}.dReadoutFOV twix_obj.hdr.MeasYaps.sSliceArray.asSlice{1}.dPhaseFOV twix_obj.hdr.MeasYaps.sSliceArray.asSlice{1}.dThickness].';
FOVxyz = FOVrps(permutedims);

hostVoxDim_mm = FOVxyz./Hxyz;
hostVoxDim_mm_rps = hostVoxDim_mm(permutedims_toXYZ);

kspaceCentre_rps = [ceil(twix_obj.image.centerCol(1)/2) twix_obj.image.centerLin(1) twix_obj.image.centerPar(1)];


nc = twix_obj.image.NCha;
nS = twix_obj.image.NSet; % for MP2RAGE this will be 2...


% and put stuff into html report and the screen
if nS == 2
    fprintf(fidHTML,['<h4>Host MP2RAGE ']);
    fprintf(['\n\n\nHost MP2RAGE ']);
else
    fprintf(fidHTML,['<h4>Host MPRAGE ']);
    fprintf(['\n\n\nHost MPRAGE ']);
end
fprintf(fidHTML,['sequence:</h4>\n']);
fprintf(['sequence:\n']);

fprintf(fidHTML,['Scanner field strength: ' num2str(twix_obj.hdr.MeasYaps.sTXSPEC.asNucleusInfo{1,1}.lFrequency/42.58e6,'%.2f') ' T<br><br>\n']);

fprintf(fidHTML,['RPS dimensions in raw data: ' num2str(hrps(1)) 'x' num2str(hrps(2)) 'x' num2str(hrps(3)) '<br>\n']);
fprintf(['RPS dimensions in raw data: ' num2str(hrps(1)) 'x' num2str(hrps(2)) 'x' num2str(hrps(3)) '\n']);
fprintf(fidHTML,['XYZ dimensions reconstructed: ' num2str(Hxyz(1)) 'x' num2str(Hxyz(2)) 'x' num2str(Hxyz(3)) '<br>\n']);
fprintf(['XYZ dimensions reconstructed: ' num2str(Hxyz(1)) 'x' num2str(Hxyz(2)) 'x' num2str(Hxyz(3)) '\n']);
fprintf(fidHTML,['FOV - ' num2str(FOVxyz(1),'%.1f') 'x' num2str(FOVxyz(2),'%.1f') 'x' num2str(FOVxyz(3),'%.1f') 'mm<br>\n']);
fprintf(['FOV - ' num2str(FOVxyz(1),'%.1f') 'x' num2str(FOVxyz(2),'%.1f') 'x' num2str(FOVxyz(3),'%.1f') 'mm\n']);
fprintf(fidHTML,['Resolution: ' num2str(hostVoxDim_mm(1),'%.3f') 'x' num2str(hostVoxDim_mm(2),'%.3f') 'x' num2str(hostVoxDim_mm(3),'%.3f') 'mm<br>\n']);
fprintf(['Resolution: ' num2str(hostVoxDim_mm(1),'%.3f') 'x' num2str(hostVoxDim_mm(2),'%.3f') 'x' num2str(hostVoxDim_mm(3),'%.3f') 'mm\n']);
fprintf(fidHTML,['Position offset (dSag dCor dTra): '  num2str([rotAndShift.sPosition.dSag rotAndShift.sPosition.dCor rotAndShift.sPosition.dTra]) '<br>\n']);
fprintf(['Position offset (dSag dCor dTra): '  num2str([rotAndShift.sPosition.dSag rotAndShift.sPosition.dCor rotAndShift.sPosition.dTra]) '\n']);
fprintf(fidHTML,['Orientation normals (dSag dCor dTra): ' num2str([rotAndShift.sNormal.dSag rotAndShift.sNormal.dCor rotAndShift.sNormal.dTra]) '<br>\n']);
fprintf(['Orientation normals (dSag dCor dTra): ' num2str([rotAndShift.sNormal.dSag rotAndShift.sNormal.dCor rotAndShift.sNormal.dTra]) '\n']);
fprintf(fidHTML,['In plane rotation: ' num2str(rotAndShift.dInPlaneRot) '<br>\n']);
fprintf(['In plane rotation: ' num2str(rotAndShift.dInPlaneRot) '\n']);
fprintf(fidHTML,['Detected orientation: ' rotAndShift.baseOrientation '<br>\n']);
fprintf(['Detected orientation: ' rotAndShift.baseOrientation '\n']);
if isfield(twix_obj.hdr.MeasYaps,'sCoilSelectMeas')
    fprintf(['Coil used: ' char(twix_obj.hdr.MeasYaps.sCoilSelectMeas.aRxCoilSelectData{1}.asList{1}.sCoilElementID.tCoilID) ', with ' num2str(nc) ' channels active\n']);
    fprintf(fidHTML,['Coil used: ' char(twix_obj.hdr.MeasYaps.sCoilSelectMeas.aRxCoilSelectData{1}.asList{1}.sCoilElementID.tCoilID) ', with ' num2str(nc) ' channels active<br>\n']);
else
    fprintf(['Coil used: ' char(twix_obj.hdr.MeasYaps.asCoilSelectMeas{1}.asList{1}.sCoilElementID.tCoilID) ', with ' num2str(nc) ' channels active\n']);
    fprintf(fidHTML,['Coil used: ' char(twix_obj.hdr.MeasYaps.asCoilSelectMeas{1}.asList{1}.sCoilElementID.tCoilID) ', with ' num2str(nc) ' channels active<br>\n']);
end
fprintf(fidHTML,'<br>\n');
fprintf('\n');


fprintf(fidHTML,'<h4>FatNav Parameters and Recon options:</h4>\n');
fprintf('FatNav Parameters and Recon options:\n');

if manualFatNavRes
    fprintf(fidHTML,['Manually selected FatNav resolution: ' num2str(reconPars.FatNavRes_mm) ' mm<br>\n']);
    fprintf(['Manually selected FatNav resolution: ' num2str(reconPars.FatNavRes_mm) ' mm\n']);
else
    fprintf(fidHTML,['Assumed FatNav resolution (based on field strength): ' num2str(reconPars.FatNavRes_mm) ' mm<br>\n']);
    fprintf(['Assumed FatNav resolution (based on field strength): ' num2str(reconPars.FatNavRes_mm) ' mm\n']);
end
fprintf(fidHTML,['bSwapFatNavHandedness: ' num2str(reconPars.bSwapFatNavHandedness) '<br>\n']);
fprintf(fidHTML,['bApplyFatNavNoseCircshift: ' num2str(reconPars.bApplyFatNavNoseCircshift) '<br>\n']);

fprintf(fidHTML,['flipAxes_xyz (aesthetic effect only) : [' num2str(reconPars.flipAxes_xyz) ']<br>\n']);
fprintf(['flipAxes_xyz (aesthetic effect only) :' num2str(reconPars.flipAxes_xyz)  '\n']);
fprintf(fidHTML,['extraFlipMat1 (it seems that this should be diag([-1 1 -1]) in most cases):<br>&nbsp;&nbsp;&nbsp;&nbsp;' num2str(reconPars.extraFlipMat1(1,:))...
                               '<br>&nbsp;&nbsp;&nbsp;&nbsp;' num2str(reconPars.extraFlipMat1(2,:))...
                               '<br>&nbsp;&nbsp;&nbsp;&nbsp;' num2str(reconPars.extraFlipMat1(3,:)) '<br>\n']);
if any(reshape(eye(4),16,1)~=reconPars.extraFlipMat2(:))
   fprintf(fidHTML,['extraFlipMat2:<br>&nbsp;&nbsp;&nbsp;&nbsp;' num2str(reconPars.extraFlipMat2(1,:))...
                                  '<br>&nbsp;&nbsp;&nbsp;&nbsp;' num2str(reconPars.extraFlipMat2(2,:))...
                                  '<br>&nbsp;&nbsp;&nbsp;&nbsp;' num2str(reconPars.extraFlipMat2(3,:))...
                                  '<br>&nbsp;&nbsp;&nbsp;&nbsp;' num2str(reconPars.extraFlipMat2(4,:)) '<br>\n']);
end
fprintf(fidHTML,['extraPositionOffsetSignFlips (it seems that this should be [1 -1 1] in most cases) : [' num2str(reconPars.extraPositionOffsetSignFlips) ']<br>\n']);

%% Check if using HEADNECK_64 receive coil, and discard channels over the neck if this is the case (would be nice to know what Siemens does...)
%%% It seems that this manual selection was only valid for one acquisition
%%% (or scanner?!) and doesn't seem to generalise to all 64-channel
%%% datasets, so functionality is removed as default behaviour. Oct 2019
% if isfield(twix_obj.hdr.MeasYaps,'sCoilSelectMeas') ...
%         && strcmp(twix_obj.hdr.MeasYaps.sCoilSelectMeas.aRxCoilSelectData{1}.asList{1}.sCoilElementID.tCoilID,'"HeadNeck_64"')
%     
%     iC_keep = 1:nc;
%     if nc==52
%         iC_keep([1 7 8 18 29 30 39 40 49 50]) = []; % these channels cover the neck (in one test dataset with 52 data channels from the 64-channel coil...)
%         fprintf(['\n\n****  Detected use of HeadNeck_64 RF coil **** \n'...
%             'Using manually predefined set of channels\n'...
%             'to reduce signal from neck area\n' ...
%             '*********************************\n\n']);
%         fprintf(fid,['<br><br>\n\n****  Detected use of HeadNeck_64 RF coil **** <br>\n'...
%             'Using manually predefined set of channels<br>\n'...
%             'to reduce signal from neck area<br>\n' ...
%             '*********************************<br><br>\n\n']);
%     end
% else
    iC_keep = 1:nc;
% end

nc_keep = length(iC_keep);
    




%% Check the number of FatNavs available compared to the size of the host data

% find the right (chronological) order for the k-space lines in the raw data

if reconPars.bLinParSwap % MP2RAGE sequence currently doesn't allow 2D GRAPPA, so this means the number of FatNavs should already match that dimension
    alignDim_rps = 3;
    iSamp = 1:hrps(alignDim_rps);           
else
    alignDim_rps = 2;
    thisLine = squeeze(twix_obj.imageWithRefscan(1,1,:,1,1,1,1,1,1,1,1,1));
    iSamp = find(thisLine);
end

if isfield(twix_obj,'FatNav')
    nFatNavs = twix_obj.FatNav.dataSize(9)/twix_obj.image.NAve;
else
    nFatNavs = 0; % presumably wasn't actually a FatNav scan!
end

fprintf(fidHTML,['<br><br><br>No. of FatNavs: ' num2str(nFatNavs) '<br>\n']);
fprintf(fidHTML,['No. of measured lines in host sequence: ' num2str(length(iSamp)) '<br>\n']);
fprintf(['\n\nNo. of FatNavs: ' num2str(nFatNavs) '\n']);
fprintf(['No. of measured lines in host sequence: ' num2str(length(iSamp)) '\n']);

if length(iSamp)~=nFatNavs
    disp('Error: number of FatNavs doesn''t seem to match number of acquired lines found')
end

%% Process the FatNavs
% - First reconstruct each FatNav, then co-register using SPM to obtain
%   motion-estimates

if nFatNavs > 0
    [ACSims, timingReport_FatNavs, fatnavdir] = processFatNavs_GRAPPA4x4(twix_obj, ...
        outDir,'FatNavRes_mm',reconPars.FatNavRes_mm, 'iAve', reconPars.iAve, ...
        'appendString', appendString,'bSwapHandedness',reconPars.bSwapFatNavHandedness,...
        'bApplyNoseCircshift',reconPars.bApplyFatNavNoseCircshift,...
        'bUseNeckMasking',reconPars.bApplyFatNavNeckMasking);
    
    % And put stuff into the html report
    if exist(fatnavdir,'dir') % could have just kept the motion-parameters file...
        imdims = 2*[304 128];
        fprintf(fidHTML,['<h4>FatNavs</h4>\n']);
        
        if exist([fatnavdir '/a_FatNav_MoCoPars_' MIDstr '.png'],'file')
            copyfile([fatnavdir '/a_FatNav_MoCoPars_' MIDstr '.png'],[htmlDir '/motion_parameters.png']);
            fprintf(fidHTML,['Estimated motion parameters:<br>\n']);
            fprintf(fidHTML,['<img src="motion_parameters.png"><br><br>\n']);
        end
        if exist([fatnavdir '/a_FatNav_ACSim_' MIDstr '.png'],'file')
            copyfile([fatnavdir '/a_FatNav_ACSim_' MIDstr '.png'],[htmlDir '/ACSim.png']);
            fprintf(fidHTML,['<br><br>FatNav ACS image:<br>\n']);
            fprintf(fidHTML,['<img src="ACSim.png" height=%s width=%s><br><br>\n'],num2str(imdims(2)),num2str(imdims(1)));
        end
        if exist([fatnavdir '/a_FatNav_NoseCircshift_' MIDstr '.png'],'file')
            copyfile([fatnavdir '/a_FatNav_NoseCircshift_' MIDstr '.png'],...
                [htmlDir '/a_FatNav_NoseCircshift_' MIDstr '.png']);
            fprintf(fidHTML,['The FatNavs have been automatically circularly shifted in the y-direction to attempt to account for a non-isocentre head position (seemingly common at ultra-high fields) which causes the nose to wrap and could affect apparent motion parameters:<br>\n']);
            fprintf(fidHTML,['<img src="a_FatNav_NoseCircshift_' MIDstr '.png"><br><br>\n']);
        end
        if exist([fatnavdir '/a_FatNav_ACSim_' MIDstr '_channelCut.png'],'file')
            copyfile([fatnavdir '/a_FatNav_ACSim_' MIDstr '_channelCut.png'],[htmlDir '/ACSim_channelCut.png']);
            fprintf(fidHTML,['ACS image (after removing neck channels from 64ch coil):<br>\n']);
            fprintf(fidHTML,['<img src="ACSim_channelCut.png" height=%s width=%s><br><br>\n'],num2str(imdims(2)),num2str(imdims(1)));
        end
        if exist([fatnavdir '/a_mov_eachFatNav.gif'],'file')
            copyfile([fatnavdir '/a_mov_eachFatNav.gif'],[htmlDir '/mov_eachFatNav.gif'])
            fprintf(fidHTML,['15 example FatNavs covering complete scan:<br>\n']);
            fprintf(fidHTML,['<img src="mov_eachFatNav.gif" height=%s width=%s><br><br>\n'],num2str(imdims(2)),num2str(imdims(1)));
        end
        if exist([fatnavdir '/a_mov_spm_eachFatNav.gif'],'file')
            copyfile([fatnavdir '/a_mov_spm_eachFatNav.gif'],[htmlDir '/mov_spm_eachFatNav.gif'])
            fprintf(fidHTML,['And after rigid-body registration using SPM:<br>\n']);
            fprintf(fidHTML,['<img src="mov_spm_eachFatNav.gif" height=%s width=%s><br><br>\n'],num2str(imdims(2)),num2str(imdims(1)));
        end
        if exist([fatnavdir '/a_FatNav_NeckMask_' MIDstr '.png'],'file')
            copyfile([fatnavdir '/a_FatNav_NeckMask_' MIDstr '.png'],...
                [htmlDir '/a_FatNav_NeckMask_' MIDstr '.png']);
            fprintf(fidHTML,['An automatic masking was applied to attempt to restrict registration to rigid regions of the FatNav:<br>\n']);
            fprintf(fidHTML,['<img src="a_FatNav_NeckMask_' MIDstr '.png" height=%s width=%s><br><br>\n'],num2str(imdims(2)),num2str(imdims(1)*2));
        end
        
        
    end
end

%% Make iFFT of host ACS for rapid orientation check without waiting for GRAPPA to finish

if Arps(2) > 1
    kdataHostACS = squeeze(twix_obj.refscan());
    kdataHostACS = kdataHostACS(:,:,:,:,end); % keep only last Set
    kdataHostACS = permute(kdataHostACS,[1 3 4 2]); % move coils to last dim
    
    imACS_rps = ssos(ifft3s(kdataHostACS));
    
    imACS_xyz = permute(imACS_rps,[permutedims 4]); % permute rest according to permutedims

    imACS_xyz = flipAxes(imACS_xyz,reconPars.flipAxes_xyz);
    
    hf = figure('Visible','off'); % make an invisible figure for all figure plots
    % to try to avoid stealing focus
    set(hf,'Position',[    50   50   950  340])
    orthoview(imACS_rps,'useNewFig',0);
    subplot1(1); axis normal; title('Read/Phase')
    ylabel('Low-res host image in native RPS coords')
    subplot1(2); axis normal; title('Read/Slice')
    subplot1(3); axis normal; title('Phase/Slice')
    export_fig([htmlDir '/orientationCheck_Host_RPS.png']);
    clf
    orthoview(imACS_xyz,'useNewFig',0);
    subplot1(1); axis normal; 
    ylabel('Low-res host image in assumed XYZ coords')
    subplot1(2); axis normal;
    subplot1(3); axis normal; 
    export_fig([htmlDir '/orientationCheck_Host_XYZ.png']);

    close(hf);
    
    fprintf(fidHTML,['<br><br><br>Orientation check for host images RPS to XYZ:<br>\n']);
    fprintf(fidHTML,['<img src="orientationCheck_Host_XYZ.png"><br><br>\n']);
    fprintf(fidHTML,['<img src="orientationCheck_Host_RPS.png"><br><br>\n']);

end

%% load moco parameters and align their orientation to the host data

if nFatNavs > 0
    fitResult = load([outDir '/motion_parameters_spm_' MIDstr appendString '.mat']);
    
    this_fitMat = fitResult.MPos_cent.mats;
    
    % Account for the fact that the host sequence may not have been acquired at
    % isocentre:
    %%% danielg: June 26 - attempt to 'fix' this code using moveFrame approach
    A_fatnav2host = eye(4);
    % add rotations:
    A_fatnav2host(1:3,1:3) = rotAndShift.RotMat*reconPars.extraFlipMat1;
    % and translations:
    A_fatnav2host(1:3,4) = -rotAndShift.RotMat*reconPars.extraFlipMat1*...
        (reconPars.extraPositionOffsetSignFlips(:).*rotAndShift.Shifts_SagCorTra(:));
    
    % it doesn't make sense to me to create a new matrix here, but when I was
    % debugging I was getting desparate. Leaving in here 'just in case'...
    A_fatnav2host_forMats = reconPars.extraFlipMat2*A_fatnav2host;
    
    this_fitMat_mm = moveFrame(this_fitMat,A_fatnav2host_forMats);
    
    if exist([fatnavdir '/eachFatNav_001.nii'],'file')
        
        %
        fileTest1 = [fatnavdir '/test1.nii'];
        fileTest2 = [fatnavdir '/test2.nii'];
        copyfile([fatnavdir '/eachFatNav_001.nii'],fileTest1);
        copyfile([fatnavdir '/eachFatNav_001.nii'],fileTest2);
        
        % attempt to use SPM directly to put FatNav into Host RPS space
        V_source = spm_vol_nifti(fileTest1);
        V_virtual = spm_vol_nifti(fileTest2);
        
        % also use niftiread functions to get header as this is (slightly)
        % different to SPM's header...
        nii = load_untouch_nii(fileTest1);
        
        
        % Change the dimensions (number of voxels) to alter your FOV
        V_virtual.dim = round(FOVrps.' / reconPars.FatNavRes_mm);
        
        %  Calculate centering shift caused purely by the FOV dimension change
        orig_center  = V_source.dim / 2;
        new_center   = V_virtual.dim / 2;
        shift_voxels = ( new_center - orig_center);
        
        % Convert voxel shift to mm using the original voxel spacing
        shift_mm = reconPars.FatNavRes_mm * shift_voxels';
        
        % Account for any nosecircshift in the fatnav
        shift_nose_mm = (nii.hdr.hist.srow_y(4) + ((V_source.dim(2)/2)*reconPars.FatNavRes_mm)); % srow_y is -ve, so this finds diff
        %     A_shiftNose = eye(4);
        %     A_shiftNose(2,4) = -shift_nose_mm;
        
        %     V_source.mat = A_fatnav2host * A_shiftNose * V_source.mat;
        
        V_source.mat = A_fatnav2host * V_source.mat;
        
        V_source.mat(1:3, 4) = V_source.mat(1:3, 4) + shift_mm;
        
        V_source.mat(2, 4) = V_source.mat(2, 4) + shift_nose_mm;
        
        % Pass both headers to spm_reslice.
        spm_reslice([V_virtual, V_source], struct('mask',false,'mean',false,'interp',1,'which',1,'wrap',[1 1 1],'prefix',''));
        
        newFatIm = rn(fileTest1);
        
        hf = figure('Visible','off');
        
        oOut = orthoview(newFatIm,'useNewFig',0);
        %
        set(gcf,'Position',[    50   50   950  540])
        subplot(oOut.hAx(1))
        title('Read/Phase')
        ylabel('First FatNav realigned to Host RPS')
        subplot(oOut.hAx(2))
        title('Read/Slice')
        subplot(oOut.hAx(3))
        title('Phase/Slice')
        export_fig([htmlDir '/orientationCheck_FatVolume.png'])
        %
        close(hf);
        
        fprintf(fidHTML,['Orientation check for host sequence slice rotation and positioning:<br>\n']);
        fprintf(fidHTML,['<img src="orientationCheck_FatVolume.png"><br><br>\n']);
        fprintf(fidHTML,['(The fat volume shown above should approximately correspond to the FOV chosen for the host sequence in the "native RPS" coords shown above that)<br><br><br>\n']);
    end
    disp('Ims done')
    
end

%% Do GRAPPA (if necessary) for the host sequence
% Note that this can use rather a lot of hard-disk space, as the raw data
% has to be shuffled around a few times to end up with one file per 
% reconstructed RF channel


if Arps(2) > 1
    if reconPars.bGetGRAPPA_SVD
        [~,~,Vsvd] = svd(reshape(kdataHostACS,[],size(kdataHostACS,4)),'econ');
        svdCombinePars = Vsvd(:,1);
    else
        svdCombinePars = [];
    end
    useGRAPPAforHost = 1;
    if reconPars.bGRAPPAinRAM        
        [grappaRecon_1DFFT, mOutGRAPPA] = performHostGRAPPArecon_RAMonly(twix_obj,...
            struct('iAve',reconPars.iAve,'iRep',reconPars.iRep),2,svdCombinePars);
        mOutGRAPPA.grappaRecon_1DFFT = grappaRecon_1DFFT; clear grappaRecon_1DFFT;
        timingReport_hostRecon = mOutGRAPPA.timingReport;
    else
%         [mOutGRAPPA, timingReport_hostRecon] =
%         performHostGRAPPArecon(twix_obj,tempDir); % the original low-RAM
%         alternative code needs verifying as the output ends up corrupted
%         (in at least some cases)
        [timingReport_hostRecon, tempNameRoots] = performHostGRAPPArecon_toDisk(twix_obj,tempDir,struct('iAve',reconPars.iAve,'iRep',reconPars.iRep),2,svdCombinePars);    
    end
else
    useGRAPPAforHost = 0;
    timingReport_hostRecon = 0;
end




%% Double-check the orientation with host data against the images from the FatNavs
% left/right matching
% - find the coil which has the biggest asymmetry left > right, and see if it
% is on the same side in both the FatNavs and the host data

if nFatNavs > 0
    
    nc_FatNavs = size(ACSims,4);
    asymData_left = sum(reshape(abs(ACSims(1:floor(FatNavDims_xyz(1)/2),:,:,:)),[],nc_FatNavs),1);
    asymData_right = sum(reshape(abs(ACSims(ceil(FatNavDims_xyz(1)/2):FatNavDims_xyz(1),:,:,:)),[],nc_FatNavs),1);
    asymFactor = asymData_left./asymData_right;
    iAsymCoil = find(asymFactor==max(asymFactor),1);
    
    if useGRAPPAforHost
        if reconPars.bGRAPPAinRAM
            hostExampleVolume = squeeze(mOutGRAPPA.grappaRecon_1DFFT(:,iAsymCoil,:,:,nS));
        else
            hostExampleVolume = zeros(hrps.');
            for iReadSlice = 1:hrps(1) % virtual 'slices' in the readout direction
                tempData = load([tempNameRoots.grappaRecon_1DFFT '_' num2str(iReadSlice) '_1_' num2str(nS) '.mat']);
                hostExampleVolume(iReadSlice,:,:) = reshape(tempData.outData(1,iAsymCoil,:,:),[1 hrps(2) hrps(3)]);
            end
        end
        hostExampleVolume = ifft1s(ifft1s(hostExampleVolume,2),3);
    else
        hostExampleVolume = squeeze(twix_obj.image(:,iAsymCoil,:,:,1,reconPars.iAve,1,1,reconPars.iRep,nS));
        hostExampleVolume = ifft3s(hostExampleVolume);
    end
    
    
    % ov1 = orthoview(ACSims(:,:,:,iAsymCoil),'drawIms',0,'mip',1,'clims',[0 max(abs(ACSims(:)))]);
    % ov2 = orthoview(hostExampleVolume,'drawIms',0,'mip',0,'clims',[0 max(abs(hostExampleVolume(:)))]);
    
    % Use SPM technique to put FatNav ACS image (single asymmetric coil only)
    % into Host RPS space:
    fileTest1 = [fatnavdir '/test1.nii'];
    fileTest2 = [fatnavdir '/test2.nii'];
    sn(abs(ACSims(:,:,:,iAsymCoil)),fileTest1,reconPars.FatNavRes_mm*[1 1 1])
    
    copyfile(fileTest1,fileTest2);
    
    % attempt to use SPM directly to put FatNav into Host RPS space
    V_source = spm_vol_nifti(fileTest1);
    V_virtual = spm_vol_nifti(fileTest2);
    
    % Change the dimensions (number of voxels) to alter your FOV
    %     V_virtual.dim = FatNavDims_HostRPS;
    V_virtual.dim = round(FOVrps.' / reconPars.FatNavRes_mm);
    
    %  Calculate centering shift caused purely by the FOV dimension change
    orig_center  = V_source.dim / 2;
    new_center   = V_virtual.dim / 2;
    shift_voxels = ( new_center - orig_center);
    
    % Convert voxel shift to mm using the original voxel spacing
    shift_mm = reconPars.FatNavRes_mm * shift_voxels';
    
    % Note that ACSim wouldn't have the nose circshift applied, so no need to
    % correct here as above
    V_source.mat = A_fatnav2host * V_source.mat;
    V_source.mat(1:3, 4) = V_source.mat(1:3, 4) + shift_mm;
    
    % Pass both headers to spm_reslice.
    spm_reslice([V_virtual, V_source], struct('mask',false,'mean',false,'interp',1,'which',1,'wrap',[0 0 0],'prefix',''));
    
    
    newFatImAsym = rn(fileTest1);
    
    ov1 = orthoview(newFatImAsym,'drawIms',0,'mip',1,'clims',[0 max(abs(ACSims(:)))]);
    ov2 = orthoview(hostExampleVolume,'drawIms',0,'mip',1,'clims',[0 max(abs(hostExampleVolume(:)))]);
    
    hf = figure('Visible','off'); % make an invisible figure for all figure plots
    % to try to avoid stealing focus
    set(hf,'Position',[    22   594   702   473])
    hAx = subplot1(2,1,'figHandle',hf);
    subplot(hAx(1))
    imab(ov1.oneIm,[0 .5*max(ov1.oneIm(:))])
    ylabel({'MIP of FatNav', 'in Host RPS space' , ['coil ' num2str(iAsymCoil)]})
    subplot(hAx(2))
    imab(ov2.oneIm,[0 .5*max(ov2.oneIm(:))])
    ylabel({'MIP of Host', 'GRAPPA recon' , [ 'coil ' num2str(iAsymCoil)]})
    export_fig([htmlDir '/orientationCheck_xy.png']);
    close(hf);
    
    fprintf(fidHTML,['Orientation check for left/right symmetry:<br>\n']);
    fprintf(fidHTML,['<img src="orientationCheck_xy.png"><br><br>\n']);
    fprintf(fidHTML,['(Both images above should have the brightest signal in the same places on the images. If not, the orientation of the FatNavs is not correctly aligned with the host sequence)<br><br><br>\n']);
    
    
    
    %% FatNavs are not acquired concurrently with the MP2RAGE data, so in the
    % case of 'brittle' motion it may be especially beneficial to average
    % neighbouring motion estimates rather than taking the one from the same TR
    % (which, in the pulse sequence, is actually closer in time to the
    % following MP2RAGE data)
    
    rotTrans = rotmat2euler(this_fitMat_mm(1:3,1:3,:));
    rotTrans(4:6,:) = squeeze(this_fitMat_mm(1:3,4,:));
    
    newRotTrans = [rotTrans(:,1) (rotTrans(:,1:end-1)+rotTrans(:,2:end))/2];
    
    this_fitMat_mm = zeros(size(this_fitMat_mm));
    this_fitMat_mm(4,4,:) = 1;
    
    this_fitMat_mm(1:3,1:3,:) = euler2rmat(newRotTrans(1:3,:));
    this_fitMat_mm(1:3,4,:) = newRotTrans(4:6,:);
    
    
    
    %% If GRAPPA was used in the 'slow' PE direction then the motion parameters will need to be interpolated to have values throughout k-space
    
    
    
    % if accelerated in fatnav direction, then motion needs to be interpolated
    % for the gaps
    if nFatNavs ~= hrps(alignDim_rps)
        rotTrans = rotmat2euler(this_fitMat_mm(1:3,1:3,:));
        rotTrans(4:6,:) = squeeze(this_fitMat_mm(1:3,4,:));
        
        rotTrans_interp = interp1(iSamp,rotTrans.',1:hrps(alignDim_rps)).';
        startIndex = iSamp(1);
        % copy in values at edges of k-space to avoid extrapolation errors
        if startIndex > 1
            rotTrans_interp(:,1:startIndex-1) = repmat(rotTrans_interp(:,startIndex),[1 startIndex-1]);
        end
        endIndex = iSamp(end);
        if endIndex < hrps(alignDim_rps)
            rotTrans_interp(:,endIndex+1:end) = repmat(rotTrans_interp(:,endIndex),[1 hrps(alignDim_rps)-endIndex]);
        end
        
        this_fitMat_mm_interp = zeros(4,4,hrps(alignDim_rps));
        this_fitMat_mm_interp(4,4,:) = 1;
        
        this_fitMat_mm_interp(1:3,1:3,:) = euler2rmat(rotTrans_interp(1:3,:));
        this_fitMat_mm_interp(1:3,4,:) = rotTrans_interp(4:6,:);
        
        fitMats_mm_toApply = this_fitMat_mm_interp;
        
    else
        fitMats_mm_toApply = this_fitMat_mm;
        fitMats_mm_toApply(4,4,:) = 1;
    end
    
    
    
    % make it so that centre of k-space is not 'moved' (accounting for partial Fourier):
    fitMats_mm_toApply = recentre_affmats(fitMats_mm_toApply,kspaceCentre_rps(alignDim_rps));
    alignIndices = 1:hrps(alignDim_rps);
    
end

%% Prepare for the retrospective motion-correction of the host sequence
tStart_applyMoco = clock;

timingReport_applyMoco = zeros(nc,nS); % store the time to reconstruct each volume with the NUFFT

fprintf('............\n')
if nFatNavs > 0
    fprintf('... Performing NUFFT on each volume to apply motion correction parameters\n');
end


%%% Declare all the variables

if ~reconPars.bFullParforRecon 
    
    if ~reconPars.bKeepReconInRAM
        tempFile = [tempDir '/tempReconData_' MIDstr '.mat'];
        if exist(tempFile,'file')
            delete(tempFile)
        end
        mOut = matfile(tempFile,'writable',true);
    end
    
    if nS==2
        mOut.all_uniImage = zeros(Hxyz(1),Hxyz(2),Hxyz(3),'single');
        mOut.all_uniImage_corrected = zeros(Hxyz(1),Hxyz(2),Hxyz(3),'single');
        mOut.all_refImage = zeros(Hxyz(1),Hxyz(2),Hxyz(3),'single');
        mOut.all_refImage_corrected = zeros(Hxyz(1),Hxyz(2),Hxyz(3),'single');
    end
    
    if strcmp(reconPars.coilCombineMethod,'lowres')        
        mOut.all_ims_ref = zeros(Hxyz(1),Hxyz(2),Hxyz(3),1,'single'); % ref is only based on INV2 as it has more SNR
        mOut.all_ims_ref_corrected = zeros(Hxyz(1),Hxyz(2),Hxyz(3),1,'single');
        
        mOut.all_ims = complex(zeros(Hxyz(1),Hxyz(2),Hxyz(3),nS,'single'));
        mOut.all_ims_corrected = complex(zeros(Hxyz(1),Hxyz(2),Hxyz(3),nS,'single'));
        
    else
        mOut.all_ims = zeros(Hxyz(1),Hxyz(2),Hxyz(3),nS,'single');
        mOut.all_ims_corrected = zeros(Hxyz(1),Hxyz(2),Hxyz(3),nS,'single');
    end
    
else % for Parfor to work you can't use the mOut structure trick (which could either be in memory or a matfile...)
    
    all_uniImage = zeros(Hxyz(1),Hxyz(2),Hxyz(3),'single');
    all_uniImage_corrected = zeros(Hxyz(1),Hxyz(2),Hxyz(3),'single');
    all_refImage = zeros(Hxyz(1),Hxyz(2),Hxyz(3),'single');
    all_refImage_corrected = zeros(Hxyz(1),Hxyz(2),Hxyz(3),'single');
    
    if strcmp(reconPars.coilCombineMethod,'lowres')
        all_ims_ref = zeros(Hxyz(1),Hxyz(2),Hxyz(3),1,'single'); % ref is only based on INV2 as it has more SNR
        all_ims_ref_corrected = zeros(Hxyz(1),Hxyz(2),Hxyz(3),1,'single');
        
        all_ims = complex(zeros(Hxyz(1),Hxyz(2),Hxyz(3),nS,'single'));
        all_ims_corrected = complex(zeros(Hxyz(1),Hxyz(2),Hxyz(3),nS,'single'));
    else
        all_ims = zeros(Hxyz(1),Hxyz(2),Hxyz(3),nS,'single');
        all_ims_corrected = zeros(Hxyz(1),Hxyz(2),Hxyz(3),nS,'single');
        
        all_ims_ref = []; % parfor seems to sniff around in unreachable parts of code and get annoyed if unused variables don't exist yet...
        all_ims_ref_corrected = [];
    end
    
    if useGRAPPAforHost % separate the two sets of MP2RAGE to allow parfor
        data_set1 = mOutGRAPPA.grappaRecon_1DFFT(:,iC_keep,:,:,1);
        if nS>1
            data_set2 = mOutGRAPPA.grappaRecon_1DFFT(:,iC_keep,:,:,2);
        else
            data_set2 = zeros(1,nc_keep); % try to keep later parfor's happy if there is no set 2
        end
    else
        data_set1 = squeeze(twix_obj.image(:,iC_keep,:,:,1,reconPars.iAve,1,1,reconPars.iRep,1));
        if nS>1
            data_set2 = squeeze(twix_obj.image(:,iC_keep,:,:,1,reconPars.iAve,1,1,reconPars.iRep,2));
        else
            data_set2 = zeros(1,nc_keep); % try to keep later parfor's happy if there is no set 2
        end
    end
    
    
end

%% Apply the motion-correction

if ~reconPars.bFullParforRecon 

    if nFatNavs > 0
        % generate the st object for NUFFT and precalculate the phase offsets
        if reconPars.bUseGPU
            [~,FT,~, phaseTranslations] = applyRetroMC_gpunufft(zeros(hrps'),fitMats_mm_toApply,alignDim_rps,alignIndices,11,...
                hostVoxDim_mm_rps,Hrps,kspaceCentre_rps,-1,reconPars.NUFFTosf,1);
        else
            [~,st,~, phaseTranslations] = applyRetroMC_nufft(zeros(hrps'),fitMats_mm_toApply,alignDim_rps,alignIndices,11,...
                hostVoxDim_mm_rps,Hrps,kspaceCentre_rps,-1,reconPars.NUFFTosf,1);
        end
    end
    
    for iC = 1:nc_keep           
        
        mOut.thiscoil_ims = complex(zeros(Hxyz(1),Hxyz(2),Hxyz(3),nS,'single'));
        mOut.thiscoil_ims_corrected = complex(zeros(Hxyz(1),Hxyz(2),Hxyz(3),nS,'single'));
        
        for iS = 1:nS
        
            tic
            fprintf(['Reconstructing coil ' num2str(iC) ' of ' num2str(nc_keep) ', set ' num2str(iS) '\n']);
            
            if useGRAPPAforHost
                if reconPars.bGRAPPAinRAM 
                    thisData = squeeze(mOutGRAPPA.grappaRecon_1DFFT(:,iC_keep(iC),:,:,iS));
                else
                    thisData = zeros(hrps.');
                    for iReadSlice = 1:hrps(1) % virtual 'slices' in the readout direction
                        tempData = load([tempNameRoots.grappaRecon_1DFFT '_' num2str(iReadSlice) '_1_' num2str(iS) '.mat']);
                        thisData(iReadSlice,:,:) = reshape(tempData.outData(1,iC_keep(iC),:,:),[1 hrps(2) hrps(3)]);
                    end
                end
                
                % use this line to have more brain coverage for debugging:
                % thisData = squeeze(mOutGRAPPA.dataCombined(:,:,:,2)); iC=1;iS = 1; 
                
                thisData = fft1s(thisData,1); % put into full 3D k-space
            else
                thisData = squeeze(twix_obj.image(:,iC_keep(iC),:,:,1,reconPars.iAve,1,1,reconPars.iRep,iS));
            end
                                
            newData = thisData;
      
%             thisData_corrected = applyRetroMC_nufft(newData,fitMats_mm_toApply,alignDim,alignIndices,11,hostVoxDim_mm,Hxyz,kspaceCentre_xyz);
%             % the multi-CPU (parfor) version of the NUFFT can be called with a
%             % negative value for the 'useTable' input:
%             %    thisData_corrected = applyRetroMC_nufft(newData,fitMats_mm_toApply,alignDim,alignIndices,-11,hostVoxDim_mm,Hxyz,kspaceCentre_xyz);
%             % However, in my latest tests (September 2016) I found that it was
%             % actually slower for the current NUFFT parameters than the
%             % non-parallelized version, so I've disabled it again as a
%             % default...
            
            if nFatNavs > 0
                if reconPars.bUseGPU
                    thisData_corrected = FT' * (newData(:).*phaseTranslations(:));
                    disp('GPU version')
                else
                    thisData_corrected =  nufft_adj_single(newData(:).*phaseTranslations(:),st);
                end
            end
            
            % apply transforms to corrected and uncorrected data for final matrix shape:
            thisData = permute(thisData,permutedims);
            if nFatNavs > 0
                thisData_corrected = permute(thisData_corrected,permutedims);
            end
                        
            if any(hxyz~=Hxyz)
                thisData(Hxyz(1),Hxyz(2),Hxyz(3)) = 0; % extend to new size
                thisData = circshift(thisData,double([Hxyz(1)-hxyz(1) Hxyz(2)-hxyz(2) Hxyz(3)-hxyz(3)]));
                % the 'double' in the above line appears necessary in certain
                % versions of Matlab, no idea why...
            end                   

            thisData = ifft3s(thisData)*prod(hxyz);

            % important to apply the flip after the FFT for full agreement
            % between iFFT and NUFFT pipelines
            thisData = flipAxes(thisData,reconPars.flipAxes_xyz);
            if nFatNavs > 0
                thisData_corrected = flipAxes(thisData_corrected,reconPars.flipAxes_xyz);
            end
            
%             save('svdCombineData.mat','thisData','newData','thisData_corrected','Hxyz','hxyz','Hrps','hrps',...
%                 'fitMats_mm_toApply','alignDim_rps','alignIndices','hostVoxDim_mm_rps','kspaceCentre_rps',...
%                 'reconPars','permutedims','-v7.3');
            
%             addpath ~/matlab/generaltools/
%             SliceBrowser2(cat(4,nma(abs(thisData)),nma(abs(thisData_corrected))))
% SliceBrowser2(cat(4,nma(abs(test1)),nma(abs(thisData_corrected))))

            
            timingReport_applyMoco(iC,iS) = toc;
            
            if nS > 1
                mOut.thiscoil_ims(:,:,:,iS) = thisData;
                if nFatNavs > 0                    
                    mOut.thiscoil_ims_corrected(:,:,:,iS) = thisData_corrected;
                end
            else
                mOut.thiscoil_ims = thisData;
                if nFatNavs > 0
                    mOut.thiscoil_ims_corrected = thisData_corrected;
                end
            end
            
        end
        
        % create coil-combined INV1 (and INV2 if available)
        switch reconPars.coilCombineMethod
            case 'default'
                mOut.all_ims = mOut.all_ims + abs(mOut.thiscoil_ims).^2;
                mOut.all_ims_corrected = mOut.all_ims_corrected + abs(mOut.thiscoil_ims_corrected).^2;
            case 'lowres'
                if nS==2
                    weightsImage = mOut.thiscoil_ims(:,:,:,2);
                    weightsImage = ifft3s(tukeyfilt3d(fft3s(weightsImage),0,filterFrac,1))*prod(hxyz);
                    mOut.all_ims = mOut.all_ims + bsxfun(@times,mOut.thiscoil_ims,conj(weightsImage));
                    mOut.all_ims_ref = mOut.all_ims_ref + abs(weightsImage).^2;
                    
                    weightsImage = mOut.thiscoil_ims_corrected(:,:,:,2);
                    weightsImage = ifft3s(tukeyfilt3d(fft3s(weightsImage),0,filterFrac,1))*prod(hxyz);
                    mOut.all_ims_corrected = mOut.all_ims_corrected + bsxfun(@times,mOut.thiscoil_ims_corrected,conj(weightsImage));
                    mOut.all_ims_ref_corrected = mOut.all_ims_ref_corrected + abs(weightsImage).^2;
                else % nS==1
                    weightsImage = mOut.thiscoil_ims;
                    weightsImage = ifft3s(tukeyfilt3d(fft3s(weightsImage),0,filterFrac,1))*prod(hxyz);
                    mOut.all_ims = mOut.all_ims + mOut.thiscoil_ims.*conj(weightsImage);
                    mOut.all_ims_ref = mOut.all_ims_ref + abs(weightsImage).^2;
                    
                    weightsImage = mOut.thiscoil_ims_corrected;
                    weightsImage = ifft3s(tukeyfilt3d(fft3s(weightsImage),0,filterFrac,1))*prod(hxyz);
                    mOut.all_ims_corrected = mOut.all_ims_corrected + mOut.thiscoil_ims_corrected.*conj(weightsImage);
                    mOut.all_ims_ref_corrected = mOut.all_ims_ref_corrected + abs(weightsImage).^2;
                end
                
                
                
        end
        
        
        % create UNI image for MP2RAGE
        if nS==2
            
            switch reconPars.coilCombineMethod
                case 'default'
                    mOut.all_uniImage = mOut.all_uniImage + ...
                        real(mOut.thiscoil_ims(:,:,:,2).*conj(mOut.thiscoil_ims(:,:,:,1))) .* (abs(mOut.thiscoil_ims(:,:,:,2)).^2) ./ ...
                        (abs(mOut.thiscoil_ims(:,:,:,1)).^2 + abs(mOut.thiscoil_ims(:,:,:,2).^2));
                    mOut.all_refImage = mOut.all_refImage + abs(mOut.thiscoil_ims(:,:,:,2)).^2;
                    
                    mOut.all_uniImage_corrected = mOut.all_uniImage_corrected + ...
                        real(mOut.thiscoil_ims_corrected(:,:,:,2).*conj(mOut.thiscoil_ims_corrected(:,:,:,1))).*(abs(mOut.thiscoil_ims_corrected(:,:,:,2)).^2) ./ ...
                        (abs(mOut.thiscoil_ims_corrected(:,:,:,1)).^2 + abs(mOut.thiscoil_ims_corrected(:,:,:,2).^2));
                    mOut.all_refImage_corrected = mOut.all_refImage_corrected + abs(mOut.thiscoil_ims_corrected(:,:,:,2)).^2;
                    
                case 'lowres'
                    weightsImage = mOut.thiscoil_ims(:,:,:,2);
                    weightsImage = ifft3s(tukeyfilt3d(fft3s(weightsImage),0,filterFrac,1))*prod(hxyz);
                    mOut.all_uniImage = mOut.all_uniImage + ...
                        real(mOut.thiscoil_ims(:,:,:,2).*conj(mOut.thiscoil_ims(:,:,:,1))) .* (abs(weightsImage).^2) ./ ...
                        (abs(mOut.thiscoil_ims(:,:,:,1)).^2 + abs(mOut.thiscoil_ims(:,:,:,2).^2));
                    mOut.all_refImage = mOut.all_refImage + abs(weightsImage).^2;
                    
                    weightsImage = mOut.thiscoil_ims_corrected(:,:,:,2);
                    weightsImage = ifft3s(tukeyfilt3d(fft3s(weightsImage),0,filterFrac,1))*prod(hxyz);
                    mOut.all_uniImage_corrected = mOut.all_uniImage_corrected + ...
                        real(mOut.thiscoil_ims_corrected(:,:,:,2).*conj(mOut.thiscoil_ims_corrected(:,:,:,1))).*(abs(weightsImage).^2) ./ ...
                        (abs(mOut.thiscoil_ims_corrected(:,:,:,1)).^2 + abs(mOut.thiscoil_ims_corrected(:,:,:,2).^2));
                    mOut.all_refImage_corrected = mOut.all_refImage_corrected + abs(weightsImage).^2;
                    
            end
            
        end
        
        % allow saving of motion-corrected complex data
        if reconPars.bKeepComplexImageData
            save([outDir '/complexImageData_coil' num2str(iC,'%.2d') '.mat'],'mOut','-v7.3');
        end
          
        
    end
    
    switch reconPars.coilCombineMethod
        case 'default'
            mOut.all_ims = sqrt(mOut.all_ims);
            mOut.all_ims_corrected = sqrt(mOut.all_ims_corrected);
        case 'lowres'
            if nS==2
                mOut.all_ims = bsxfun(@rdivide,abs(mOut.all_ims),sqrt(mOut.all_ims_ref));
                mOut.all_ims_corrected = bsxfun(@rdivide,abs(mOut.all_ims_corrected),sqrt(mOut.all_ims_ref_corrected));
            else
                mOut.all_ims = abs(mOut.all_ims)./sqrt(mOut.all_ims_ref);
                mOut.all_ims_corrected = abs(mOut.all_ims_corrected)./sqrt(mOut.all_ims_ref_corrected);
            end
            
    end
    
    if nS==2
        mOut.all_uniImage = mOut.all_uniImage./mOut.all_refImage;
        sn(int16(4095*mOut.all_ims(:,:,:,1)/max(reshape(mOut.all_ims,[],1))),[outDir '/INV1'],hostVoxDim_mm)
        sn(int16(4095*(mOut.all_uniImage+0.5)) ,[outDir '/UNI'],hostVoxDim_mm)
        sn(int16(4095*mOut.all_ims(:,:,:,2)/max(reshape(mOut.all_ims,[],1))),[outDir '/INV2'],hostVoxDim_mm)
        
        if nFatNavs > 0
            mOut.all_uniImage_corrected = mOut.all_uniImage_corrected./mOut.all_refImage_corrected;
            sn(int16(4095*mOut.all_ims_corrected(:,:,:,1)/max(reshape(mOut.all_ims_corrected,[],1))),[outDir '/INV1_' outputSuffix],hostVoxDim_mm)
            sn(int16(4095*(mOut.all_uniImage_corrected+0.5)),[outDir '/UNI_' outputSuffix],hostVoxDim_mm)
            sn(int16(4095*mOut.all_ims_corrected(:,:,:,2)/max(reshape(mOut.all_ims_corrected,[],1))),[outDir '/INV2_' outputSuffix],hostVoxDim_mm)
        end
        
       
    else
        % if using matfile variables you can't use index in dimensions which
        % aren't there, even if you only put a 1 there...!
        sn(int16(4095*mOut.all_ims/max(reshape(mOut.all_ims,[],1))),[outDir '/INV1'],hostVoxDim_mm)
        if nFatNavs > 0
            sn(int16(4095*mOut.all_ims_corrected/max(reshape(mOut.all_ims_corrected,[],1))),[outDir '/INV1_' outputSuffix],hostVoxDim_mm)
        end
        
    end
    
    
else % the much faster version with much hungrier RAM requirements:
    
    if nFatNavs > 0
    % generate the st object for NUFFT and precalculate the phase offsets    
    [~,st,~, phaseTranslations] = applyRetroMC_nufft(zeros(hrps'),fitMats_mm_toApply,alignDim_rps,alignIndices,11,...
        hostVoxDim_mm_rps,Hrps,kspaceCentre_rps,-1,reconPars.NUFFTosf,1);
    else
        phaseTranslations = []; % parfor is weird and it seems these variables must exist even if not used
        st = [];
    end
    
    parfor iC = 1:nc_keep
        
        thiscoil_ims = complex(zeros(Hxyz(1),Hxyz(2),Hxyz(3),nS,'single'));
        thiscoil_ims_corrected = complex(zeros(Hxyz(1),Hxyz(2),Hxyz(3),nS,'single'));

        
        for iS = 1:nS
        
            fprintf(['Reconstructing coil ' num2str(iC) ' of ' num2str(nc_keep) ', set ' num2str(iS) '\n']);
                        
            switch iS
                case 1
                    thisData = squeeze(data_set1(:,iC,:,:));
                case 2
                    thisData = squeeze(data_set2(:,iC,:,:));
            end
            if useGRAPPAforHost
                thisData = fft1s(thisData,1); % put into full 3D k-space
            end
            
                
            tic       
            if nFatNavs > 0
                thisData_corrected =  nufft_adj_single(thisData(:).*phaseTranslations(:),st);
            end
            timingReport_applyMoco(iC,iS) = toc;
            
            % apply transforms to corrected and uncorrected data for final matrix shape:
            thisData = permute(thisData,permutedims);
            if nFatNavs > 0
                thisData_corrected = permute(thisData_corrected,permutedims);
            end
            
            if any(hxyz~=Hxyz)
                thisData(Hxyz(1),Hxyz(2),Hxyz(3)) = 0; % extend to new size
                thisData = circshift(thisData,double([Hxyz(1)-hxyz(1) Hxyz(2)-hxyz(2) Hxyz(3)-hxyz(3)]));
                % the 'double' in the above line appears necessary in certain
                % versions of Matlab, no idea why...
            end
            
            thisData = ifft3s(thisData)*prod(hxyz);
            
            % important to apply the flip after the FFT for full agreement
            % between iFFT and NUFFT pipelines
            thisData = flipAxes(thisData,reconPars.flipAxes_xyz);
            if nFatNavs > 0
                thisData_corrected = flipAxes(thisData_corrected,reconPars.flipAxes_xyz);
            end
            
            if nS > 1
                thiscoil_ims(:,:,:,iS) = thisData;
                if nFatNavs > 0
                    thiscoil_ims_corrected(:,:,:,iS) = thisData_corrected;
                end
            else
                thiscoil_ims = thisData;
                if nFatNavs > 0
                    thiscoil_ims_corrected = thisData_corrected;
                end
            end
            
        end
        
        % create coil-combined INV1 (and INV2 if available)
        switch reconPars.coilCombineMethod
            case 'default'
                all_ims = all_ims + abs(thiscoil_ims).^2;
                all_ims_corrected = all_ims_corrected + abs(thiscoil_ims_corrected).^2;
            case 'lowres'
                if nS==2
                    weightsImage = thiscoil_ims(:,:,:,2);
                    weightsImage = ifft3s(tukeyfilt3d(fft3s(weightsImage),0,filterFrac,1))*prod(hxyz);
                    all_ims = all_ims + bsxfun(@times,thiscoil_ims,conj(weightsImage));
                    all_ims_ref = all_ims_ref + abs(weightsImage).^2;
                    
                    weightsImage = thiscoil_ims_corrected(:,:,:,2);
                    weightsImage = ifft3s(tukeyfilt3d(fft3s(weightsImage),0,filterFrac,1))*prod(hxyz);
                    all_ims_corrected = all_ims_corrected + bsxfun(@times,thiscoil_ims_corrected,conj(weightsImage));
                    all_ims_ref_corrected = all_ims_ref_corrected + abs(weightsImage).^2;
                else % nS==1
                    weightsImage = thiscoil_ims;
                    weightsImage = ifft3s(tukeyfilt3d(fft3s(weightsImage),0,filterFrac,1))*prod(hxyz);
                    all_ims = all_ims + thiscoil_ims.*conj(weightsImage);
                    all_ims_ref = all_ims_ref + abs(weightsImage).^2;
                    
                    weightsImage = thiscoil_ims_corrected;
                    weightsImage = ifft3s(tukeyfilt3d(fft3s(weightsImage),0,filterFrac,1))*prod(hxyz);
                    all_ims_corrected = all_ims_corrected + thiscoil_ims_corrected.*conj(weightsImage);
                    all_ims_ref_corrected = all_ims_ref_corrected + abs(weightsImage).^2;
                end
                
        end
        
        
        % create UNI image for MP2RAGE
        if nS==2
            
            switch reconPars.coilCombineMethod
                case 'default'
                    all_uniImage = all_uniImage + ...
                        real(thiscoil_ims(:,:,:,2).*conj(thiscoil_ims(:,:,:,1))) .* (abs(thiscoil_ims(:,:,:,2)).^2) ./ ...
                        (abs(thiscoil_ims(:,:,:,1)).^2 + abs(thiscoil_ims(:,:,:,2).^2));
                    all_refImage = all_refImage + abs(thiscoil_ims(:,:,:,2)).^2;
                    
                    all_uniImage_corrected = all_uniImage_corrected + ...
                        real(thiscoil_ims_corrected(:,:,:,2).*conj(thiscoil_ims_corrected(:,:,:,1))).*(abs(thiscoil_ims_corrected(:,:,:,2)).^2) ./ ...
                        (abs(thiscoil_ims_corrected(:,:,:,1)).^2 + abs(thiscoil_ims_corrected(:,:,:,2).^2));
                    all_refImage_corrected = all_refImage_corrected + abs(thiscoil_ims_corrected(:,:,:,2)).^2;
                    
                case 'lowres'
                    weightsImage = thiscoil_ims(:,:,:,2);
                    weightsImage = ifft3s(tukeyfilt3d(fft3s(weightsImage),0,filterFrac,1))*prod(hxyz);
                    all_uniImage = all_uniImage + ...
                        real(thiscoil_ims(:,:,:,2).*conj(thiscoil_ims(:,:,:,1))) .* (abs(weightsImage).^2) ./ ...
                        (abs(thiscoil_ims(:,:,:,1)).^2 + abs(thiscoil_ims(:,:,:,2).^2));
                    all_refImage = all_refImage + abs(weightsImage).^2;
                    
                    weightsImage = thiscoil_ims_corrected(:,:,:,2);
                    weightsImage = ifft3s(tukeyfilt3d(fft3s(weightsImage),0,filterFrac,1))*prod(hxyz);
                    all_uniImage_corrected = all_uniImage_corrected + ...
                        real(thiscoil_ims_corrected(:,:,:,2).*conj(thiscoil_ims_corrected(:,:,:,1))).*(abs(weightsImage).^2) ./ ...
                        (abs(thiscoil_ims_corrected(:,:,:,1)).^2 + abs(thiscoil_ims_corrected(:,:,:,2).^2));
                    all_refImage_corrected = all_refImage_corrected + abs(weightsImage).^2;
            end
            
        end
        
        % allow saving of motion-corrected complex data
        if reconPars.bKeepComplexImageData
            parsave([outDir '/complexImageData_coil' num2str(iC,'%.2d') '.mat'],thiscoil_ims,thiscoil_ims_corrected);
        end
        
        
    end % parfor
    
    
    switch reconPars.coilCombineMethod
        case 'default'
            all_ims = sqrt(all_ims);
            all_ims_corrected = sqrt(all_ims_corrected);
        case 'lowres'
            if nS==2
                all_ims = bsxfun(@rdivide,abs(all_ims),sqrt(all_ims_ref));
                all_ims_corrected = bsxfun(@rdivide,abs(all_ims_corrected),sqrt(all_ims_ref_corrected));
            else
                all_ims = abs(all_ims)./sqrt(all_ims_ref);
                all_ims_corrected = abs(all_ims_corrected)./sqrt(all_ims_ref_corrected);
            end
    end
    
    
    all_uniImage = all_uniImage./all_refImage;
    
    sn(int16(4095*all_ims(:,:,:,1)/max(reshape(all_ims,[],1))),[outDir '/INV1'],hostVoxDim_mm)
    
    if nFatNavs > 0
        all_uniImage_corrected = all_uniImage_corrected./all_refImage_corrected;
        sn(int16(4095*all_ims_corrected(:,:,:,1)/max(reshape(all_ims_corrected,[],1))),[outDir '/INV1_' outputSuffix],hostVoxDim_mm)
    end
    
    if nS>1
        sn(int16(4095*(all_uniImage+0.5)) ,[outDir '/UNI'],hostVoxDim_mm)
        sn(int16(4095*all_ims(:,:,:,2)/max(reshape(all_ims,[],1))),[outDir '/INV2'],hostVoxDim_mm)
        if nFatNavs > 0
            sn(int16(4095*all_ims_corrected(:,:,:,2)/max(reshape(all_ims_corrected,[],1))),[outDir '/INV2_' outputSuffix],hostVoxDim_mm)
            sn(int16(4095*(all_uniImage_corrected+0.5)),[outDir '/UNI_' outputSuffix],hostVoxDim_mm)
        end
    end
    
    mOut.all_ims = all_ims;
    clear all_ims
    mOut.all_ims_corrected = all_ims_corrected;
    clear all_ims_corrected
    mOut.all_uniImage = all_uniImage;
    clear all_uniImage;
    mOut.all_uniImage_corrected = all_uniImage_corrected;
    clear all_uniImage_corrected;
    
end





%%


fprintf('Done\n');


tFinish_applyMoco = clock;
timingReport_totalTimeApplyMoco = etime(tFinish_applyMoco,tStart_applyMoco);
avgTimeApplyMocoPerVolume = sum(timingReport_applyMoco(:))/ nc / nS;

%% And put the reconstructed images into the html

fprintf(fidHTML,'<h4>Reconstructed MP(2)RAGE images:</h4>\n');

switch nS

    case 1 % again different code because matfile variables can't index dimensions which don't exist
        
        ov1 = orthoview(mOut.all_ims,'drawIms',0);
        imab_overwrite([htmlDir '/INV1.png'],ov1.oneIm);
        if nFatNavs > 0            
            ov1 = orthoview(mOut.all_ims_corrected,'drawIms',0);
            imab_overwrite([htmlDir '/INV1_' outputSuffix '.png'],ov1.oneIm);
        end
        
        fprintf(fidHTML,['INV1 image before correction:<br>\n']);
        fprintf(fidHTML,['<img src="INV1.png"><br><br>\n']);
        if nFatNavs > 0
            fprintf(fidHTML,['INV1 image after correction:<br>\n']);
            fprintf(fidHTML,['<img src="INV1_' outputSuffix '.png"><br><br>\n']);
        end
        
        if nFatNavs > 0
            
            testMagick = system('convert -version');
            
            if testMagick==0 % can use ImageMagick to make animated GIFs...
                processString = ['convert -dispose 2 -delay 50 -loop 0 ' htmlDir '/INV1*.png ' htmlDir '/mov_INV1.gif'];
                system(processString);
                fprintf(fidHTML,['INV1 image movie before/after correction:<br>\n']);
                fprintf(fidHTML,['<img src="mov_INV1.gif"><br><br>\n']);
            else
                testMagickNew = system('magick -version'); % test for newer version of imageMagick (thanks to Sila Dokumaci!)
                if testMagickNew==0
                    processString = ['magick -dispose 2 -delay 50 -loop 0 ' htmlDir '/INV1*.png ' htmlDir '/mov_INV1.gif'];
                    system(processString);
                    fprintf(fidHTML,['INV1 image movie before/after correction:<br>\n']);
                    fprintf(fidHTML,['<img src="mov_INV1.gif"><br><br>\n']);
                end
                
            end
        end
        
        %%% Add a sagittal zoom of the images before and after correction
        xi = round(Hxyz(1)/2+Hxyz(1)/10);
        yi = round(Hxyz(2)/2+Hxyz(2)/8: .9*Hxyz(2)); % arbitrary cut off points for zoom in y and z
        zi = round(Hxyz(3)/2+Hxyz(3)/8: .8*Hxyz(3));

        clim1 = percentile(mOut.all_ims,97);
        clim1c = percentile(mOut.all_ims_corrected,97);
        zoomLimsScale = 1.3; % arbitrary factor as zoomed images tend to be clipped
        
        hf = figure('Visible','off');
        set(hf,'Position',[   246   611   500   500])
        imab(squeeze(mOut.all_ims(xi,yi,zi)),[0 clim1*zoomLimsScale])
        colormap(gray)
        export_fig([htmlDir '/zoom.png'])
        close(hf);
        
        if nFatNavs > 0
            hf = figure('Visible','off');
            set(hf,'Position',[   246   611   500   500])
            imab(squeeze(mOut.all_ims_corrected(xi,yi,zi)),[0 clim1c*zoomLimsScale])
            colormap(gray)
            export_fig([htmlDir '/zoom_' outputSuffix '.png'])
            close(hf);
        end
        
        if nFatNavs > 0
            
            if testMagick==0
                processString = ['convert -dispose 2 -delay 50 -loop 0 ' htmlDir '/zoom.png ' htmlDir '/zoom_' outputSuffix '.png ' htmlDir '/mov_zoom.gif'];
                system(processString);
                
                fprintf(fidHTML,['Sagittal zoom before/after correction:<br>\n']);
                fprintf(fidHTML,['<img src="mov_zoom.gif"><br><br>\n']);
            elseif testMagickNew==0
                processString = ['magick -dispose 2 -delay 50 -loop 0 ' htmlDir '/zoom.png ' htmlDir '/zoom_' outputSuffix '.png ' htmlDir '/mov_zoom.gif'];
                system(processString);
                
                fprintf(fidHTML,['Sagittal zoom before/after correction:<br>\n']);
                fprintf(fidHTML,['<img src="mov_zoom.gif"><br><br>\n']);
            else
                fprintf(fidHTML,['Sagittal zoom before correction:<br>\n']);
                fprintf(fidHTML,['<img src="zoom.png"><br><br>\n']);
                fprintf(fidHTML,['Sagittal after correction:<br>\n']);
                fprintf(fidHTML,['<img src="zoom_' outputSuffix '.png"><br><br>\n']);
            end
        end
        
    case 2
        
        ov1 = orthoview(mOut.all_ims(:,:,:,1),'drawIms',0);
        imab_overwrite([htmlDir '/INV1.png'],ov1.oneIm);
        if nFatNavs > 0            
            ov1 = orthoview(mOut.all_ims_corrected(:,:,:,1),'drawIms',0);
            imab_overwrite([htmlDir '/INV1_' outputSuffix '.png'],ov1.oneIm);
        end
        
        fprintf(fidHTML,['INV1 image before correction:<br>\n']);
        fprintf(fidHTML,['<img src="INV1.png"><br><br>\n']);
        if nFatNavs > 0
            fprintf(fidHTML,['INV1 image after correction:<br>\n']);
            fprintf(fidHTML,['<img src="INV1_' outputSuffix '.png"><br><br>\n']);
        end
        
        if nFatNavs > 0
            
            testMagick = system('convert -version');
            
            if testMagick==0 % can use ImageMagick to make animated GIFs...
                processString = ['convert -dispose 2 -delay 50 -loop 0 ' htmlDir '/INV1*.png ' htmlDir '/mov_INV1.gif'];
                system(processString);
                fprintf(fidHTML,['INV1 image movie before/after correction:<br>\n']);
                fprintf(fidHTML,['<img src="mov_INV1.gif"><br><br>\n']);
                testMagickNew = -1; % create this variable in case it needs to exist for code to run
            else
                testMagickNew = system('magick -version');  % test for newer version of imageMagick (thanks to Sila Dokumaci!)
                if testMagickNew==0
                    processString = ['magick -dispose 2 -delay 50 -loop 0 ' htmlDir '/INV1*.png ' htmlDir '/mov_INV1.gif'];
                    system(processString);
                    fprintf(fidHTML,['INV1 image movie before/after correction:<br>\n']);
                    fprintf(fidHTML,['<img src="mov_INV1.gif"><br><br>\n']);
                end
                
            end
        end
        
        ov1 = orthoview(mOut.all_ims(:,:,:,2),'drawIms',0);
        imab_overwrite([htmlDir '/INV2.png'],ov1.oneIm);

        ov1 = orthoview(mOut.all_uniImage,'drawIms',0);
        imab_overwrite([htmlDir '/UNI.png'],ov1.oneIm);
        
        if nFatNavs > 0
            ov1 = orthoview(mOut.all_ims_corrected(:,:,:,2),'drawIms',0);
            imab_overwrite([htmlDir '/INV2_' outputSuffix '.png'],ov1.oneIm);
            ov1 = orthoview(mOut.all_uniImage_corrected,'drawIms',0);
            imab_overwrite([htmlDir '/UNI_' outputSuffix '.png'],ov1.oneIm);
        end
        
        fprintf(fidHTML,['INV2 image before correction:<br>\n']);
        fprintf(fidHTML,['<img src="INV2.png"><br><br>\n']);
        if nFatNavs > 0            
            fprintf(fidHTML,['INV2 image after correction:<br>\n']);
            fprintf(fidHTML,['<img src="INV2_' outputSuffix '.png"><br><br>\n']);
        end
        
        if nFatNavs > 0
            if testMagick==0
                processString = ['convert -dispose 2 -delay 50 -loop 0 ' htmlDir '/INV2.png ' htmlDir '/INV2_' outputSuffix '.png ' htmlDir '/mov_INV2.gif'];
                system(processString);
                processString = ['convert -dispose 2 -delay 50 -loop 0 ' htmlDir '/UNI.png ' htmlDir '/UNI_' outputSuffix '.png ' htmlDir '/mov_UNI.gif'];
                system(processString);
                fprintf(fidHTML,['INV2 image movie before/after correction:<br>\n']);
                fprintf(fidHTML,['<img src="mov_INV2.gif"><br><br>\n']);
                
            elseif testMagickNew==0
                processString = ['magick -dispose 2 -delay 50 -loop 0 ' htmlDir '/INV2.png ' htmlDir '/INV2_' outputSuffix '.png ' htmlDir '/mov_INV2.gif'];
                system(processString);
                processString = ['magick -dispose 2 -delay 50 -loop 0 ' htmlDir '/UNI.png ' htmlDir '/UNI_' outputSuffix '.png ' htmlDir '/mov_UNI.gif'];
                system(processString);
                fprintf(fidHTML,['INV2 image movie before/after correction:<br>\n']);
                fprintf(fidHTML,['<img src="mov_INV2.gif"><br><br>\n']);
            end
        end
        
        fprintf(fidHTML,['UNI image before correction:<br>\n']);
        fprintf(fidHTML,['<img src="UNI.png"><br><br>\n']);
        if nFatNavs > 0
            fprintf(fidHTML,['UNI image after correction:<br>\n']);
            fprintf(fidHTML,['<img src="UNI_' outputSuffix '.png"><br><br>\n']);
        end
        
        if nFatNavs > 0
            if testMagick==0 || testMagickNew==0
                fprintf(fidHTML,['UNI image movie before/after correction:<br>\n']);
                fprintf(fidHTML,['<img src="mov_UNI.gif"><br><br>\n']);
            end
        end
        
        
        if nFatNavs > 0
            
            % Do skull stripping with BET to be able to see the vessels within the brain more
            % clearly in the INV2 image
            disp('...............')
            disp('... Checking for FSL...')
            testFSL = getenv('FSLDIR');
            if ~isempty(testFSL)
                disp('... Found FSL, assuming that ''bet'' and ''fslmaths'' commands will work')
                disp('... Performing BET brain extraction ...')
                setenv('FSLOUTPUTTYPE','NIFTI'); % this uses more space than NIFTI_GZ, but stays compatible with SPM
                system(['bet ' outDir '/INV2_' outputSuffix ' ' outDir '/INV2_' outputSuffix '_bet -m -f 0.2']);
                system(['fslmaths ' outDir '/INV2 -mul ' outDir '/INV2_' outputSuffix '_bet_mask ' outDir '/INV2_bet']);
                disp('... Done')
                
                %%% make MIPs of INV2 after BET
                
                inv2 = rn([outDir '/INV2_bet.nii']);
                inv2c = rn([outDir '/INV2_' outputSuffix '_bet.nii']);
                
                ov1 = orthoview(inv2,'mip',1,'drawIms',0);
                ov2 = orthoview(inv2c,'mip',1,'drawIms',0);
                
                im1 = abs(ov1.oneIm); im1 = im1/max(im1(:));
                im2 = abs(ov2.oneIm); im2 = im2/max(im2(:));
                
                imab_overwrite([htmlDir '/INV2_MIP.png'],im1);
                imab_overwrite([htmlDir '/INV2_' outputSuffix '_MIP.png'],im2);
                
                if testMagick==0
                    processString = ['convert -dispose 2 -delay 50 -loop 0 ' htmlDir '/INV2_MIP.png ' htmlDir '/INV2_' outputSuffix '_MIP.png ' htmlDir '/mov_INV2_MIP.gif'];
                    system(processString);
                    
                    fprintf(fidHTML,['INV2 MIP movie before/after correction:<br>\n']);
                    fprintf(fidHTML,['<img src="mov_INV2_MIP.gif"><br><br>\n']);
                elseif testMagickNew==0
                    processString = ['magick -dispose 2 -delay 50 -loop 0 ' htmlDir '/INV2_MIP.png ' htmlDir '/INV2_' outputSuffix '_MIP.png ' htmlDir '/mov_INV2_MIP.gif'];
                    system(processString);
                    
                    fprintf(fidHTML,['INV2 MIP movie before/after correction:<br>\n']);
                    fprintf(fidHTML,['<img src="mov_INV2_MIP.gif"><br><br>\n']);
                else
                    fprintf(fidHTML,['INV2 MIP before correction:<br>\n']);
                    fprintf(fidHTML,['<img src="INV2_MIP.png"><br><br>\n']);
                    fprintf(fidHTML,['INV2 MIP after correction:<br>\n']);
                    fprintf(fidHTML,['<img src="INV2_' outputSuffix '_MIP.png"><br><br>\n']);
                end
                
            end
            
            %%% Add a zoom of the 3 images before and after correction
            if ~isempty(testFSL) % there will be a BET mask to use for choosing the brain region
                brainmask = rn([outDir '/INV2_' outputSuffix '_bet_mask.nii']);
                xi = round(Hxyz(1)/2+Hxyz(1)/10);
                yi = round(Hxyz(2)/2+Hxyz(2)/8:find(squeeze(any(any(brainmask,1),3)),1,'last'));
                zi = round(Hxyz(3)/2:find(squeeze(any(any(brainmask,1),2)),1,'last'));
            else
                xi = round(Hxyz(1)/2+Hxyz(1)/10);
                yi = round(Hxyz(2)/2+Hxyz(2)/8: .9*Hxyz(2)); % arbitrary cut off points for zoom in y and z
                zi = round(Hxyz(3)/2+Hxyz(3)/8: .8*Hxyz(3));
            end
            clim1 = percentile(mOut.all_ims(:,:,:,1),97);
            clim2 = percentile(mOut.all_ims(:,:,:,2),97);
            clim1c = percentile(mOut.all_ims_corrected(:,:,:,1),97);
            clim2c = percentile(mOut.all_ims_corrected(:,:,:,2),97);
            clims_uni = [-.5 .5];
            zoomLimsScale = 1.3; % arbitrary factor as zoomed images tend to be clipped
            
            hf = figure('Visible','off');
            set(hf,'Position',[   246   611   982   494])
            hAx = subplot1(1,3,'figHandle',hf);
            subplot(hAx(1))
            imab(squeeze(mOut.all_ims(xi,yi,zi,1)),[0 clim1*zoomLimsScale])
            subplot(hAx(2))
            imab(squeeze(mOut.all_ims(xi,yi,zi,2)),[0 clim2*zoomLimsScale])
            subplot(hAx(3))
            imab(squeeze(mOut.all_uniImage(xi,yi,zi)),clims_uni)
            colormap(gray)
            export_fig([htmlDir '/zoom.png'])
            close(hf);
            
            hf = figure('Visible','off');
            set(hf,'Position',[   246   611   982   494])
            hAx = subplot1(1,3,'figHandle',hf);
            subplot(hAx(1))
            imab(squeeze(mOut.all_ims_corrected(xi,yi,zi,1)),[0 clim1c*zoomLimsScale])
            subplot(hAx(2))
            imab(squeeze(mOut.all_ims_corrected(xi,yi,zi,2)),[0 clim2c*zoomLimsScale])
            subplot(hAx(3))
            imab(squeeze(mOut.all_uniImage_corrected(xi,yi,zi)),clims_uni)
            colormap(gray)
            export_fig([htmlDir '/zoom_' outputSuffix '.png'])
            close(hf);
            
            if testMagick==0
                processString = ['convert -dispose 2 -delay 50 -loop 0 ' htmlDir '/zoom.png ' htmlDir '/zoom_' outputSuffix '.png ' htmlDir '/mov_zoom.gif'];
                system(processString);
                
                fprintf(fidHTML,['Sagittal zoom before/after correction:<br>\n']);
                fprintf(fidHTML,['<img src="mov_zoom.gif"><br><br>\n']);
            elseif testMagickNew==0
                processString = ['magick -dispose 2 -delay 50 -loop 0 ' htmlDir '/zoom.png ' htmlDir '/zoom_' outputSuffix '.png ' htmlDir '/mov_zoom.gif'];
                system(processString);
                
                fprintf(fidHTML,['Sagittal zoom before/after correction:<br>\n']);
                fprintf(fidHTML,['<img src="mov_zoom.gif"><br><br>\n']);
            else
                fprintf(fidHTML,['Sagittal zoom before correction:<br>\n']);
                fprintf(fidHTML,['<img src="zoom.png"><br><br>\n']);
                fprintf(fidHTML,['Sagittal after correction:<br>\n']);
                fprintf(fidHTML,['<img src="zoom_' outputSuffix '.png"><br><br>\n']);
            end
            
            
        end
        
end




%% And put a timing report into the html

stopTime = clock;
timingReport_postProcessing = etime(stopTime,tFinish_applyMoco);
totalTime = etime(stopTime,startTime)/60/60;
totalTime_hrs = floor(totalTime);
if totalTime_hrs > 0
    totalTime_mins = round(rem(totalTime,totalTime_hrs)*60);
else
    totalTime_mins = round(totalTime*60);
end

fprintf(fidHTML,['<h4>Total reconstruction time: ' num2str(totalTime_hrs) ' hours, ' num2str(totalTime_mins) ' mins</h4>\n']);
if nFatNavs > 0
    fprintf(fidHTML,['<strong>Calculate GRAPPA weights for FatNavs: </strong>' num2str(round(timingReport_FatNavs.calculateGRAPPAweights)) ' seconds.<br>\n']);
    fprintf(fidHTML,['<strong>Reconstruct FatNavs: </strong>' num2str(nFatNavs) 'x ' num2str(round(mean(timingReport_FatNavs.eachFatNav))) ' seconds. Total time (possibly parallelized!): = ' num2str(round(timingReport_FatNavs.allFatNavs)) ' seconds. <br>\n']);
    fprintf(fidHTML,['<strong>Align FatNavs using SPM: </strong>' num2str(round(timingReport_FatNavs.SPMalignment)) ' seconds.<br>\n']);
end

if reconPars.bGRAPPAinRAM
    fprintf(fidHTML,'<em>GRAPPA recon of host performed in RAM (i.e. faster)...</em><br>\n');
else
    fprintf(fidHTML,'<em>GRAPPA recon of host stored as temporary files on hard drive (i.e. slower)...</em><br>\n');
end
if useGRAPPAforHost
    fprintf(fidHTML,['<strong>Do 1D FFT for each ''slice'' in readout direction of host data:</strong> ' num2str(round(timingReport_hostRecon.FFTperSlice)) ' seconds.<br>\n']);
    fprintf(fidHTML,['<strong>Declare variables for GRAPPA recon:</strong> ' num2str(round(timingReport_hostRecon.declareVariables)) ' seconds.<br>\n']);
    fprintf(fidHTML,['<strong>GRAPPA recon:</strong> ' num2str(timingReport_hostRecon.GRAPPArecon/60,'%.1f') ' mins.<br>\n']);
end
fprintf(fidHTML,['<strong>Application of retrospective motion-correction: </strong>' num2str(nc_keep) ' channels, ' num2str(nS) ' sets, each taking ' num2str(round(avgTimeApplyMocoPerVolume)) ' seconds (possibly parallelized) = ' num2str(timingReport_totalTimeApplyMoco/60,'%.1f') ' mins.<br>\n']);



fprintf(fidHTML,'</body></html>\n');
fclose(fidHTML);


%% Delete the temporary files (which could be rather large...!)
 
% clear mOut and mOutGRAPPA files

    
if reconPars.bKeepGRAPPArecon
    if reconPars.bGRAPPAinRAM
        save([outDir '/GRAPPArecons_beforeMoco.mat'],'mOutGRAPPA','-v7.3');
    end
else
    if ~reconPars.bKeepReconInRAM
        delete(mOut.Properties.Source)
    end    
    if ~reconPars.bGRAPPAinRAM
        if ~isempty(tempNameRoots.grappaRecon_1DFFT)
            delete([tempNameRoots.grappaRecon_1DFFT '*.mat']);
        end
        if ~isempty(tempNameRoots.reconSoS)
            delete([tempNameRoots.reconSoS '*.mat']);
        end
        if ~isempty(tempNameRoots.dataCombined)
            delete([tempNameRoots.dataCombined '*.mat']);
        end
    end
    rmdir(tempDir)
end

if reconPars.bZipNIFTIs
    gzip([outDir '/*.nii']);
    delete([outDir '/*.nii']);
end

if ~reconPars.bKeepFatNavs && nFatNavs > 0
    rmdir(fatnavdir,'s')
end
    
timingReport.totalTime_hrs = totalTime_hrs;
timingReport.totalTime_mins = totalTime_mins;
if nFatNavs > 0
    timingReport.timingReport_FatNavs = timingReport_FatNavs;
    timingReport.avgTimeApplyMocoPerVolume = avgTimeApplyMocoPerVolume;
    timingReport.timingReport_totalTimeApplyMoco = timingReport_totalTimeApplyMoco;    
else
    timingReport.timingReport_FatNavs = 0;
    timingReport.avgTimeApplyMocoPerVolume = 0;
    timingReport.timingReport_totalTimeApplyMoco = 0;
end
timingReport.timingReport_hostRecon = timingReport_hostRecon;
timingReport.postProcessing = timingReport_postProcessing;

end

