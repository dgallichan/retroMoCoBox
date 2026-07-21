function timingReport = reconstructSiemensSPACEvolume(twix_obj,reconPars)
% function timingReport = reconstructSiemensSPACEvolume(twix_obj,reconPars)
% 
%
% -- July 2026, gallichand@cardiff.ac.uk -> major overhaul for v1.0.0dev
%     see retroMocoBoxVersion.m for more info

retroMocoBoxVersion = reconPars.retroMocoBoxVersion; % put this into the HTML for reference

%% Turn on 'ignoreSeg' flags for SPACE

twix_obj.image.flagIgnoreSeg = 1;
twix_obj.imageWithRefscan.flagIgnoreSeg = 1;
twix_obj.refscan.flagIgnoreSeg = 1;

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


% intialize html index file
fidHTML = fopen([htmlDir '/index.html'],'w');
fprintf(fidHTML,['<html><head><title>SPACE with FatNavs - Summary</title>\n']);
fprintf(fidHTML,'</head>\n');
fprintf(fidHTML,'<body>\n');
fprintf(fidHTML,['<h2>SPACE with FatNavs - Summary: %s</h2>\n'],[fileName]);

% include version number
fprintf(fidHTML,['<em>' char(datetime) '- created with reconstructSiemensSPACEwithFatNavs.m, version: ' retroMocoBoxVersion '</em><br><br>\n']);






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
fprintf(fidHTML,['<h4>Host SPACE sequence:</h4>\n']);
fprintf(['\n\n\nHost SPACE sequence:\n']);


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
%     iC_keep = 1:nc;
% end

% nc_keep = length(iC_keep);
    
%% check if CAIPI is used or not:

if Arps(2)>1 || Arps(3) > 1
    bIsAccelerated = 1;
    
    dataExampleMask = squeeze(twix_obj.image(1,1,:,:,1));
    dataExampleMask(dataExampleMask~=0) = 1;
    
    % First get a point close to the centre
    [xx,yy] = find(dataExampleMask);
    rr2 = (xx-kspaceCentre_rps(2)).^2 + (yy-kspaceCentre_rps(3)).^2;
    [~,iMin] = min(rr2);
    xx = xx(iMin); yy = yy(iMin);
    
    % Now check if the point offset by acceleration is sampled or not
    if dataExampleMask(xx+Arps(2),yy+Arps(3))
        bUseCAIPI = 0;
        fprintf(['Detected: GRAPPA acceleration ' num2str(Arps(2)) 'x' num2str(Arps(3)) '\n'])        ;
        fprintf(fidHTML,['Detected: GRAPPA acceleration ' num2str(Arps(2)) 'x' num2str(Arps(3)) '<br>\n']);
    else
        bUseCAIPI = 1;
        fprintf(['Detected: CAIPI acceleration ' num2str(Arps(2)) 'x' num2str(Arps(3)) '\n'])        ;
        fprintf(fidHTML,['Detected: CAIPI acceleration ' num2str(Arps(2)) 'x' num2str(Arps(3)) '<br>\n']);

    end
    

else
    bIsAccelerated = 0;
    bUseCAIPI = 0;
    fprintf(['No GRAPPA acceleration detected\n']);
    fprintf(fidHTML,['No GRAPPA acceleration detected<br>\n']);
end


%% Check the Acq order of the host sequence (SPACE does funny k-space trajectories)
% and check number of FatNavs available compared to the size of the host data

[idxMap_acrossEchoTrains, idxMap_withinEchoTrain, lineCounterMap, timeMap, rawCorrectMap] = getKspaceAcqOrder(twix_obj);

nFatNavs = twix_obj.FatNav.dataSize(9)/twix_obj.image.NAve;
nReadoutTrains = length(unique(idxMap_acrossEchoTrains(:)));


% interpolate the idxmap
[ix,iy] = find(idxMap_acrossEchoTrains);
idxMapInterp = griddata(iy,ix,idxMap_acrossEchoTrains(idxMap_acrossEchoTrains>0),(1:size(idxMap_acrossEchoTrains,2))',(1:size(idxMap_acrossEchoTrains,1)),'nearest');



hf = figure('Visible','off'); % make an invisible figure for all figure plots
% to try to avoid stealing focus
set(hf,'Position',[    50   50    1103         462])
subplot1(1,2)
subplot1(1)
imab(idxMap_acrossEchoTrains)
colorbar
title('SPACE k-space sampling, chronological')
xlabel('Phase dir')
ylabel('Slice dir')
subplot1(2)
imab(idxMapInterp)
colorbar
title('Interpolated')
xlabel('Phase dir')
ylabel('Slice dir')
export_fig([htmlDir '/SPACE_acqOrder.png']);
close(hf);

fprintf(fidHTML,['<br><br><img src="SPACE_acqOrder.png"<br><br>\n']);


fprintf(fidHTML,['<br><br><br>No. of FatNavs: ' num2str(nFatNavs) '<br>\n']);
fprintf(fidHTML,['No. of measured readout trains in host sequence: ' num2str(nReadoutTrains) '<br>\n']);
fprintf(['\n\nNo. of FatNavs: ' num2str(nFatNavs) '\n']);
fprintf(['No. of measured readout trains in host sequence: ' num2str(nReadoutTrains) '\n']);

if nReadoutTrains~=nFatNavs
    disp('Error: number of FatNavs doesn''t seem to match number of readout trains found')
end

%% Process the FatNavs
% - First reconstruct each FatNav, then co-register using SPM to obtain
%   motion-estimates


% find first non-zero fatnav in SPACEE data (why are there zeros at the start sometimes anyway...? who knows...)
testdata = squeeze(twix_obj.FatNav(1,1,1,1,1,1,1,1,:));
iFirstFatNav = find(testdata,1,'first');

[FatNav_ACSims, timingReport_FatNavs, fatnavdir] = processFatNavs_GRAPPA4x4(twix_obj, ... 
            outDir,'FatNavRes_mm',reconPars.FatNavRes_mm, 'iAve', reconPars.iAve, ...
            'appendString', appendString,'bSwapHandedness',reconPars.bSwapFatNavHandedness,...
            'bApplyNoseCircshift',reconPars.bApplyFatNavNoseCircshift,...
            'bUseNeckMasking',reconPars.bApplyFatNavNeckMasking,...
            'indexOfFirstNav',iFirstFatNav);

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


%% check orientations (for no GRAPPA this will also load all data into RAM)

if bIsAccelerated
%     Make iFFT of host ACS for rapid orientation check without waiting for GRAPPA to finish
    kdataHostACS = squeeze(twix_obj.refscan());
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

else
    
    disp('Loading data into RAM...')
    tic
    allData = twix_obj.image(); % loads everything into RAM 
    toc
    disp('Done')
    
    im_rps = squeeze(ssos(ifft1s(ifft1s(ifft1s(allData,1),3),4),2)); % coils is dim2!
   
    im_xyz = permute(im_rps,[permutedims 4]); % permute rest according to permutedims
    im_xyz = flipAxes(im_xyz,reconPars.flipAxes_xyz); 
    
    hf = figure('Visible','off'); % make an invisible figure for all figure plots
    % to try to avoid stealing focus
    set(hf,'Position',[    50   50   950  340])
    orthoview(im_rps,'useNewFig',0);
    subplot1(1); axis normal; title('Read/Phase')
    ylabel('Host image (no MoCo) in native RPS coords')
    subplot1(2); axis normal; title('Read/Slice')
    subplot1(3); axis normal; title('Phase/Slice')
    export_fig([htmlDir '/orientationCheck_Host_RPS.png']);
    clf
    orthoview(im_xyz,'useNewFig',0);
    subplot1(1); axis normal; 
    ylabel('Host image (no MoCo) in assumed XYZ coords')
    subplot1(2); axis normal;
    subplot1(3); axis normal; 
    export_fig([htmlDir '/orientationCheck_Host_XYZ.png']);

    close(hf);
    
    fprintf(fidHTML,['<br><br><br>Orientation check for host images RPS to XYZ:<br>\n']);
    fprintf(fidHTML,['<img src="orientationCheck_Host_XYZ.png"><br><br>\n']);
    fprintf(fidHTML,['<img src="orientationCheck_Host_RPS.png"><br><br>\n']);

end

%% load moco parameters and align their orientation to the host data

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

if exist([fatnavdir '/eachFatNav_' num2str(iFirstFatNav,'%.3d') '.nii'],'file')
    
   %
    fileTest1 = [fatnavdir '/test1.nii'];
    fileTest2 = [fatnavdir '/test2.nii'];
    copyfile([fatnavdir '/eachFatNav_' num2str(iFirstFatNav,'%.3d') '.nii'],fileTest1);
    copyfile([fatnavdir '/eachFatNav_' num2str(iFirstFatNav,'%.3d') '.nii'],fileTest2);
    
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
% %     A_shiftNose(2,4) = -shift_nose_mm; 
%     A_shiftNose(2,4) = shift_nose_mm;
   
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


%% Try to smooth motion parameters from measured k-space points across the points that are filled in with GRAPPA

mpars = mats2pars(this_fitMat_mm);
mpars_ims = zeros(hrps(2),hrps(3),6);
nPerSeg = ceil(twix_obj.imageWithRefscan.NAcq / twix_obj.imageWithRefscan.NSeg);

for iT = 1:nPerSeg

    [ix,iy] = find(idxMap_acrossEchoTrains==iT);
    
    for iX = 1:length(ix)
        mpars_ims(ix(iX),iy(iX),:) = reshape(mpars(:,iT),[1 1 6]);
    end
end

[ix,iy] = find(idxMap_acrossEchoTrains);
mpars_ims_smooth = zeros(hrps(2),hrps(3),6);
for iP = 1:6
    thispar = mpars_ims(:,:,iP);
    mpars_ims_smooth(:,:,iP) = griddata(iy,ix,thispar(idxMap_acrossEchoTrains>0),(1:size(idxMap_acrossEchoTrains,2))',(1:size(idxMap_acrossEchoTrains,1)),'cubic');    
end
mpars_ims_smooth(isnan(mpars_ims_smooth)) = 0;
fitMats_smoothInterp = pars2affmats(reshape(mpars_ims_smooth,hrps(2)*hrps(3),6)');
idx_forSmooth = reshape(1:hrps(2)*hrps(3),hrps(2),hrps(3));

iCentre = idxMapInterp(kspaceCentre_rps(2),kspaceCentre_rps(3)); % very difficult for accelerated TSE as 39 is adjacent to 54!
this_fitMat_mm = recentre_affmats(this_fitMat_mm,iCentre);



%% Get GRAPPA kernel for host


if bIsAccelerated    
    
    % % check if CAIPI
    % % if isfield(twix_obj.hdr.MeasYaps.sPat, 'ulCaipirinhaMode')
    % %     bUseCAIPI = twix_obj.hdr.MeasYaps.sPat.ulCaipirinhaMode;
    % % else
    % %     bUseCAIPI = 0;
    % % end
    
    % code above doesnt' work as CaipiMode can be set and used or not used it seems...
    
    
    if bUseCAIPI
        % not all kernels covered, not all code given fully tested!
        switch Arps(2)
            case 2
                if Arps(3)==2
                    % % create grapkernel for CAIPI 2x2:
                    grapKernel = zeros(5,5);
                    grapKernel([1 5],[1 3 5]) = 1;
                    grapKernel(3,[2 4]) = 1;
                    grapKernel([2 3 4],3) = 0.5;
                    
                elseif Arps(3)==3
                    % % create grapkernel for CAIPI 2x3:
                    if twix_obj.hdr.MeasYaps.sPat.lReorderingShift3D == 2
                        grapKernel = zeros(5,5);
                        grapKernel(1,[1 4]) = 1;
                        grapKernel(3,3) = 1;
                        grapKernel(5,[2 5]) = 1;
                        grapKernel(2,[2 3 4]) = 0.5;
                        grapKernel(3,[2 4]) = 0.5;
                    else
                        grapKernel = zeros(5,7);
                        grapKernel([1 3 5],[1 7]) = 1;
                        grapKernel([2 4],4) = 1;
                        grapKernel(3,[3 4 5]) = 0.5;
                        grapKernel = grapKernel';
                    end
                else
                    error('CAIPI acceleration combo not yet implemented....')
                end
                
            case 3
                if Arps(3)==3
                    %% create grapkernel for CAIPI 3x3:
                    grapKernel = zeros(6,7);
                    grapKernel([1 4],1) = 1;
                    grapKernel([2 5],4) = 1;
                    grapKernel([3 6],7) = 1;
                    grapKernel([3 4],[3 4 5]) = 0.5;
                    grapKernel(2,[3 5]) = 0.5;
                    grapKernel = grapKernel';
                    
                else
                    error('CAIPI acceleration combo not yet implemented....')
                end
                
            otherwise
                error('CAIPI acceleration combo not yet implemented....')
                % % create grapkernels for CAIPI 3x2:
                
                
                
                
        end
    else
        if Arps(2)==2 && Arps(3)==2
            grapKernel = zeros(7,7);
            grapKernel([1:2:7],[1:2:7]) = 1;
            grapKernel(4,[3 4]) = 0.5;
            grapKernel(3,4) = 0.5;
        elseif Arps(2)==2 && Arps(3)==3
            grapKernel = zeros(5,7);
            grapKernel([1:2:5],[1:3:7]) = 1;
            grapKernel(2,[3 4 5]) = 0.5;
            grapKernel(3,[3 5]) = 0.5;
        else
            error('Acceleration combo not recognised')
        end
    end
    
end

        
        
        
    
    

%% Double-check the orientation with host data against the images from the FatNavs
% left/right matching
% - find the coil which has the biggest asymmetry left > right, and see if it
% is on the same side in both the FatNavs and the host data
nc_FatNavs = size(FatNav_ACSims,4);
asymData_left = sum(reshape(abs(FatNav_ACSims(1:floor(FatNavDims_xyz(1)/2),:,:,:)),[],nc_FatNavs),1);
asymData_right = sum(reshape(abs(FatNav_ACSims(ceil(FatNavDims_xyz(1)/2):FatNavDims_xyz(1),:,:,:)),[],nc_FatNavs),1);
asymFactor = asymData_left./asymData_right;
iAsymCoil = find(asymFactor==max(asymFactor),1);

if Arps(2) > 1 || Arps(3) > 1   
    hostExampleVolume = ifft3s(kdataHostACS(:,:,:,iAsymCoil));
else
    hostExampleVolume = ifft3s(squeeze(twix_obj.image(:,iAsymCoil,:,:)));
end

% Use SPM technique to put FatNav ACS image (single asymmetric coil only)
% into Host RPS space:
fileTest1 = [fatnavdir '/test1.nii'];
fileTest2 = [fatnavdir '/test2.nii'];   
sn(abs(FatNav_ACSims(:,:,:,iAsymCoil)),fileTest1,reconPars.FatNavRes_mm*[1 1 1])

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

ov1 = orthoview(newFatImAsym,'drawIms',0,'mip',1);
ov2 = orthoview(hostExampleVolume,'drawIms',0,'mip',1);

hf = figure('Visible','off'); % make an invisible figure for all figure plots
% to try to avoid stealing focus
set(hf,'Position',[    22   594   702   473])
hAx = subplot1(2,1,'figHandle',hf);
subplot(hAx(1))
imab(ov1.oneIm,[0 .5*max(abs(FatNav_ACSims(:)))])
ylabel({'MIP of FatNav', 'in Host RPS space' , ['coil ' num2str(iAsymCoil)]})
subplot(hAx(2))
imagesc(ov2.oneIm.',[0 .5*max(abs(hostExampleVolume(:)))]); axis tight; set(gca,'xtick',[],'ytick',[])
ylabel({'MIP of Host', 'GRAPPA recon' , [ 'coil ' num2str(iAsymCoil)]})
export_fig([htmlDir '/orientationCheck_xy.png']);
close(hf);

fprintf(fidHTML,['Orientation check for left/right symmetry:<br>\n']);
fprintf(fidHTML,['<img src="orientationCheck_xy.png"><br><br>\n']);
fprintf(fidHTML,['(Both images above should have the brightest signal in the same places on the images for resolving L/R symmetry. If not, the orientation of the FatNavs is not correctly aligned with the host sequence)<br><br><br>\n']);



%% Create tukey-based edge filters to suppress strong ringing in this SPACE sequence
% this is rather arbitrary at the moment!

pFourier_PE = hrps(2)/Hrps(2);
pFourier_SL = hrps(3)/Hrps(3);



if pFourier_PE < .8
    tScale = 0.25;
    tFilt = tukey(kspaceCentre_rps(2),tScale,1-tScale);
    
    kFilt_PE = ones(hrps(2),1);
    kFilt_PE(1:round(kspaceCentre_rps(2)/2-1)) = tFilt(2:round(kspaceCentre_rps(2)/2));
else
    kFilt_PE = ones(hrps(2),1);
end

if pFourier_SL < .8
    tScale = 0.25;
    tFilt = tukey(kspaceCentre_rps(3),tScale,1-tScale);
    
    kFilt_SL = ones(1,hrps(3));
    kFilt_SL(1:round(kspaceCentre_rps(3)/2-1)) = tFilt(2:round(kspaceCentre_rps(3)/2));
else
    kFilt_SL = ones(1,hrps(3));
end

kFilt = kFilt_PE.*kFilt_SL;


if isfield(twix_obj.hdr.MeasYaps.sKSpace,'ucEnableEllipticalScanning') && twix_obj.hdr.MeasYaps.sKSpace.ucEnableEllipticalScanning
    [kx,ky] = ndgrid([1:hrps(2)] - kspaceCentre_rps(2),[1:hrps(3)] - kspaceCentre_rps(3));
    kx = kx/(hrps(2)-kspaceCentre_rps(2));
    ky = ky/(hrps(3)-kspaceCentre_rps(3));
    kr = sqrt(kx.^2+ky.^2);
    
    tFactor = 0.8;
    tFilt = tukey(hrps(3),tFactor,1-tFactor);
    tFilt = tFilt(end/2:end);
    
    krFilt = interp1(linspace(0,1,length(tFilt)),tFilt,kr);
    krFilt(isnan(krFilt)) = 0;
    
    kFilt = kFilt.*krFilt;
end




%% Do GRAPPA if necessary, otherwise just apply the tukey windowing for any partial fourier

if bIsAccelerated
    
    allData = twix_obj.image(); % loads everything into RAM (no ACS though)
    
    ACSdata = twix_obj.refscan();   
    allData_grap = zeros(size(allData));

    %% start with 1D FFT in readout direction

    allData = ifft1s(allData,1);
    ACSdata = ifft1s(ACSdata,1);
    
    %% Now do GRAPPA on each virtual 'slice' in SLC direction


    parfor iS = 1:hrps(1)
        
        % move coils to last dimension for this 'slice'
        thisData = permute(squeeze(allData(iS,:,:,:)),[2 3 1]);
        thisACS = permute(squeeze(ACSdata(iS,:,:,:)),[2 3 1]);
        
        grapW = GrappaCalib3D_arb(thisACS,grapKernel,1e-2);
        
        [iFirst_phase, iFirst_slice] = getGRAPPAstartPositions(thisData(:,:,1),grapKernel);
        
        if bUseCAIPI
            switch Arps(2)
                case 2
                    if twix_obj.hdr.MeasYaps.sPat.lReorderingShift3D == 2
                        %%% I'm not actually sure how to interpret this one, but it
                        %%% does seem to affect the 'flavour' of the CAIPI
                        %%% acceleration...
                        % triple steps in phase direction and start three times
                        grapRecon = GrappaReco3D_arb(thisData,grapKernel,grapW,...
                            [Arps(2)*3 Arps(3)],[iFirst_phase iFirst_slice;...
                            iFirst_phase+Arps(2) iFirst_slice+2; ...
                            iFirst_phase+2*Arps(2) iFirst_slice+1]);
                    else
                        % double steps in phase direction and start twice
                        grapRecon = GrappaReco3D_arb(thisData,grapKernel,grapW,...
                            [Arps(2)*2 Arps(3)],[iFirst_phase iFirst_slice;iFirst_phase+Arps(2) iFirst_slice-1+Arps(3)]);
                    end
                case 3
                    if Arps(3)==3
                        % triple steps in phase direction and start three times
                        grapRecon = GrappaReco3D_arb(thisData,grapKernel,grapW,...
                            [Arps(2)*3 Arps(3)],[iFirst_phase iFirst_slice;...
                            iFirst_phase+Arps(2) iFirst_slice+1; ...
                            iFirst_phase+2*Arps(2) iFirst_slice+2]);
                    end
            end
        else
            grapRecon = GrappaReco3D_arb(thisData,grapKernel,grapW,...
                [Arps(2) Arps(3)],[iFirst_phase iFirst_slice]);
        end
        
        
        % Put ACS back in:
        grapReconWithACS = grapRecon;
        
        % the indices here were checked manually for one particular dataset -
        % it's possible that it does not fully generalise for all scan
        % parameters!
        grapReconWithACS([1:size(thisACS,1)]+kspaceCentre_rps(2)-size(thisACS,1)/2-1,...
            [1:size(thisACS,2)]+kspaceCentre_rps(3)-size(thisACS,2)/2-1,:) = thisACS;
        
        
        %     allData_grap(iS,:,:,:) = permute(grapReconWithACS,[4 3 1 2]);
        allData_grap(iS,:,:,:) = permute(grapReconWithACS.*kFilt,[4 3 1 2]);
        
        fprintf('.');
        
    end
else
    
    allData = allData.*reshape(kFilt,[1 1 hrps(2) hrps(3)]);
    
end




 %% Create an SVD-based coil combination
% 
% [~,~,V] = svd(reshape(FatNav_ACSims,[],nc),'econ');
% combinePars = V(:,1);
% 
% %% Do coil combination to get one with good signal (almost) everywhere
% 
% grapReconSVD = zeros(hrps');
% for iC = 1:nc
%     grapReconSVD = grapReconSVD + squeeze(allData_grap(:,iC,:,:)*combinePars(iC));
% end
% 
% grapReconSVD = ifft1s(ifft1s(grapReconSVD,2),3);
% grapReconSVD = permute(grapReconSVD,permutedims);
% 
% 
% % testGrapRecon = permute(squeeze(ssos(ifft1s(ifft1s(allData_grap,3),4),2)),permutedims);
% 
% % orthoview(testGrapRecon(:,:,end:-1:1),'centre',[100 100 200])
% % orthoview(testGrapRecon(:,:,end:-1:1))


%% Do the motion correction on each coil separately


all_Image = zeros(Hxyz.');
all_Image_corrected = zeros(Hxyz.');

if bIsAccelerated
    % Undo the readout FFT to get back to k-space
    allData = fft1s(allData_grap,1);
end
    
tStart = clock;

iIdx = [2 3]; % the two phase-encoding directions

[~,st, newNewData, phaseTranslations] = applyRetroMC_nufft(zeros(hrps'),this_fitMat_mm, ...
    iIdx, idxMapInterp, 11,hostVoxDim_mm,Hrps,kspaceCentre_rps,-1,1.5,1);

tStart_applyMoco = clock;

parfor iC = 1:nc
    
    %%
    
    disp(iC)
    
    thisData = squeeze(allData(:,iC,:,:));

%     thisData = fft3s(grapReconSVD); iC = 1;% just for testing...!    
   
    newData = thisData;
    
    tic
    thisData_corrected =  nufft_adj_single(newData(:).*phaseTranslations(:),st);
    toc
    
    % apply transforms to corrected and uncorrected data for final matrix shape:
    thisData = permute(thisData,permutedims);
    thisData_corrected = permute(thisData_corrected,permutedims);
    
    if any(hxyz~=Hxyz)
        thisData(Hxyz(1),Hxyz(2),Hxyz(3)) = 0; % extend to new size
        thisData = circshift(thisData,double([Hxyz(1)-hxyz(1) Hxyz(2)-hxyz(2)-1 Hxyz(3)-hxyz(3)]));
        % the 'double' in the above line appears necessary in certain
        % versions of Matlab, no idea why...
    end
    
    thisData = ifft3s(thisData)*prod(hxyz);

    % important to apply the flip after the FFT for full agreement
    % between iFFT and NUFFT pipelines
    thisData = flipAxes(thisData,reconPars.flipAxes_xyz);
    thisData_corrected = flipAxes(thisData_corrected,reconPars.flipAxes_xyz);
    
    % Save out full complex data if desired
    if reconPars.bKeepComplexImageData
        parsave([outDir '/coil_' num2str(iC,'%.2d') '.mat'],thisData,thisData_corrected);
    end
    disp(iC)
    
    

    

    %%
    all_Image = all_Image + abs(thisData).^2;    
    all_Image_corrected = all_Image_corrected + abs(thisData_corrected).^2;
       

    
end



all_Image = sqrt(all_Image);
all_Image_corrected = sqrt(all_Image_corrected);

sn( int16(4095*(abs(all_Image)/percentile(abs(all_Image(:)),99))) ,[outDir '/SPACE'],hostVoxDim_mm)
sn( int16(4095*(abs(all_Image_corrected)/percentile(abs(all_Image_corrected(:)),99))),[outDir '/SPACE_corrected'],hostVoxDim_mm)


fprintf('Done\n');


tFinish_applyMoco = clock;
timingReport_totalTimeApplyMoco = etime(tFinish_applyMoco,tStart_applyMoco);
avgTimeApplyMocoPerVolume = sum(timingReport_totalTimeApplyMoco(:))/ nc;


%% And put the reconstructed images into the html

fprintf(fidHTML,'<h4>Reconstructed SPACE images:</h4>\n');

clims = [0 percentile(all_Image(:),97)];
clims_c = [0 percentile(all_Image_corrected(:),97)];

im1 = orthoview(all_Image,'drawIms',0);
im2 = orthoview(all_Image_corrected,'drawIms',0);

imab_overwrite([htmlDir '/SPACE_orthoview.png'],im1.oneIm,clims);
imab_overwrite([htmlDir '/SPACE_orthoview_' outputSuffix '.png'],im2.oneIm,clims_c);

fprintf(fidHTML,['SPACE image before correction:<br>\n']);
fprintf(fidHTML,['<img src="SPACE_orthoview.png"><br><br>\n']);
fprintf(fidHTML,['SPACE image after correction:<br>\n']);
fprintf(fidHTML,['<img src="SPACE_orthoview_' outputSuffix '.png"><br><br>\n']);
        
        
testMagick = system('convert -version');

if testMagick==0 % can use ImageMagick to make animated GIFs...
    processString = ['convert -dispose 2 -delay 50 -loop 0 ' htmlDir '/SPACE_orthoview.png ' htmlDir ...
                     '/SPACE_orthoview_' outputSuffix '.png ' htmlDir '/mov_SPACE.gif'];
    system(processString);
    fprintf(fidHTML,['SPACE image movie before/after correction:<br>\n']);
    fprintf(fidHTML,['<img src="mov_SPACE.gif"><br><br>\n']);
else
    testMagickNew = system('magick -version'); % test for newer version of imageMagick (thanks to Sila Dokumaci!)
    if testMagickNew==0
        processString = ['magick -dispose 2 -delay 50 -loop 0 ' htmlDir '/SPACE_orthoview.png ' htmlDir ...
                         '/SPACE_orthoview_' outputSuffix '.png ' htmlDir '/mov_SPACE.gif'];
        system(processString);
        fprintf(fidHTML,['SPACE image movie before/after correction:<br>\n']);
        fprintf(fidHTML,['<img src="mov_SPACE.gif"><br><br>\n']);
    end
    
end

%%% Add a sagittal zoom of the images before and after correction
xi = round(Hxyz(1)/2+Hxyz(1)/10);
yi = round(Hxyz(2)/2+Hxyz(2)/8: .9*Hxyz(2)); % arbitrary cut off points for zoom in y and z
zi = round(Hxyz(3)/2+Hxyz(3)/8: .8*Hxyz(3));


zoomLimsScale = 1.6; % arbitrary factor as zoomed images tend to be clipped

hf = figure('Visible','off');
set(hf,'Position',[   246   611   500   500])
imab(squeeze(all_Image(xi,yi,zi)),clims*zoomLimsScale)
colormap(gray)
export_fig([htmlDir '/zoom.png'])
close(hf);

hf = figure('Visible','off');
set(hf,'Position',[   246   611   500   500])
imab(squeeze(all_Image_corrected(xi,yi,zi)),clims_c*zoomLimsScale)
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
fprintf(fidHTML,['<strong>Calculate GRAPPA weights for FatNavs: </strong>' num2str(round(timingReport_FatNavs.calculateGRAPPAweights)) ' seconds.<br>\n']);
fprintf(fidHTML,['<strong>Reconstruct FatNavs: </strong>' num2str(nFatNavs) 'x ' num2str(round(mean(timingReport_FatNavs.eachFatNav))) ' seconds. Total time (possibly parallelized!): = ' num2str(round(timingReport_FatNavs.allFatNavs)) ' seconds. <br>\n']);
fprintf(fidHTML,['<strong>Align FatNavs using SPM: </strong>' num2str(round(timingReport_FatNavs.SPMalignment)) ' seconds.<br>\n']);


fprintf(fidHTML,'</body></html>\n');
fclose(fidHTML);



if reconPars.bZipNIFTIs
    gzip([outDir '/*.nii']);
    delete([outDir '/*.nii']);
end

if ~reconPars.bKeepFatNavs
    rmdir(fatnavdir,'s')
end
    
timingReport.totalTime_hrs = totalTime_hrs;
timingReport.totalTime_mins = totalTime_mins;
timingReport.timingReport_FatNavs = timingReport_FatNavs;
timingReport.avgTimeApplyMocoPerVolume = avgTimeApplyMocoPerVolume;
timingReport.timingReport_totalTimeApplyMoco = timingReport_totalTimeApplyMoco;
timingReport.postProcessing = timingReport_postProcessing;


fprintf('*************************************************************\n')
fprintf('***** reconstructSiemensSPACEvolume.m completed! *****\n')
fprintf('*************************************************************\n')
fprintf(['Total reconstruction time: ' num2str(totalTime_hrs) ' hours, ' num2str(totalTime_mins) ' mins\n']);


end

