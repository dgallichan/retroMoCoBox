function reconstructSiemensGREwithFatNavs_cluster(reconPars)



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

startTime = clock;


%% Make raw data object (parses file, but does not load into RAM)

tic
twix_obj = mapVBVD_fatnavs(reconPars.rawDataFile,'removeOS',1);
timingReport_parseRawDataFile = toc;

if length(twix_obj)>1 % on VE (and presumably VD as well) the raw data typically also has the adjcoilsens data as scan 1, so skip this
    twix_obj = twix_obj{2};
end

if ~isfield(twix_obj,'FatNav')
    disp('Error, no FatNavs found in raw data file!')
    return
end


%%

% for 7T-GRE sequence, the FatNav resolution is obtainable from the WIP
% parameters:
try
    reconPars.FatNavRes_mm = twix_obj.hdr.MeasYaps.sWiPMemBlock.alFree{5};
catch
    reconPars.FatNavRes_mm = 4; % if it didn't work (which it doesn't for the smaller reconcheck dataset...), assume it was 4
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


if ~isfield(reconPars,'iAve')
    reconPars.iAve = 1;
end
if ~isfield(reconPars,'iRep')
    reconPars.iRep = 1;
end

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

figIndex = 999; % figure to use for outputs

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

outDir = [reconPars.outRoot '/GRErecon_' MIDstr appendString];
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
tempDir = [reconPars.tempRoot '/temp_' getenv('SLURM_JOB_ID') '_' MIDstr appendString]; % try to include SLURM JOB ID here to avoid clashes when doing simultaneous runs
if ~exist(tempDir,'dir')
    mkdir(tempDir)
end

% intialize html index file
fidHTML = fopen([htmlDir '/index.html'],'w');
fprintf(fidHTML,['<html><head><title>GRE with FatNavs - Summary</title>\n']);
fprintf(fidHTML,'</head>\n');
fprintf(fidHTML,'<body>\n');
fprintf(fidHTML,['<h2>GRE with FatNavs - Summary: %s</h2>\n'],[fileName]);
fprintf(fidHTML,['<br>\n']);

% include version number
fprintf(fidHTML,['<em>' char(datetime) '- created with reconstructSiemensGREwithFatNavs_cluster.m, version: ' retroMocoBoxVersion '</em><br><br>\n']);





%%


% in the next line I assume that Siemens always use the same structure for
% the 'tReferenceImage0' field - but I haven't looked for documentation to
% support this, so it may not always extract the scan date properly...
if isfield(twix_obj.hdr.MeasYaps,'tReferenceImage0') % some files apparently might not even have this field...
    iDot = strfind(char(twix_obj.hdr.MeasYaps.tReferenceImage0),'.'); % some versions seem to create a char already, others not...
    thisDateString = char(twix_obj.hdr.MeasYaps.tReferenceImage0);
    if strcmp(thisDateString(iDot(end)+(1:6)),'300000') % it seems sometimes this 300000 is present, then the year is only two digits...
        fprintf(fidHTML,['<strong>Date of scan:</strong> ' local_reformatDateString(['20' thisDateString(iDot(end)+6+(1:6))]) '<br>\n']);
    else
        fprintf(fidHTML,['<strong>Date of scan:</strong> ' local_reformatDateString(thisDateString(iDot(end)+(1:8))) '<br>\n']);
    end
end
if ischar(reconPars.bKeepPatientInfo)
    fprintf(fidHTML,['<strong>Database ID:</strong> ' reconPars.bKeepPatientInfo '<br>\n']);
else
    if reconPars.bKeepPatientInfo
        %% Put patient info into HTML
        
        if isfield(twix_obj.hdr.Config,'PatientName')            
            fprintf(fidHTML,['<strong>Patient Name:</strong> ' num2str(twix_obj.hdr.Config.PatientName) '<br>\n']);
        else
            if isfield(twix_obj.hdr.Config,'tPatientName') % why would this be different for some scanners / software versions...?!?
                fprintf(fidHTML,['<strong>Patient Name:</strong> ' num2str(twix_obj.hdr.Config.tPatientName) '<br>\n']);
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
        fprintf(fidHTML,['<strong>Patient date of birth:</strong> ' local_reformatDateString(num2str(twix_obj.hdr.Config.PatientBirthDay)) '<br>\n']);
        
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
nEco = twix_obj.image.NEco;

% and put stuff into html report and the screen

fprintf(fidHTML,['<h4>Host GRE ']);
fprintf(['\n\n\nHost GRE ']);

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

if reconPars.manualFatNavRes
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
% if isfield(twix_obj.hdr.MeasYaps,'sCoilSelectMeas') ...
%         && strcmp(twix_obj.hdr.MeasYaps.sCoilSelectMeas.aRxCoilSelectData{1}.asList{1}.sCoilElementID.tCoilID,'"HeadNeck_64"')
%     iC_keep = 1:nc;
%     iC_keep([1 7 8 18 29 30 39 40 49 50]) = []; % these channels cover the neck (in one test dataset with 52 data channels from the 64-channel coil...)
%     fprintf(['\n\n****  Detected use of HeadNeck_64 RF coil **** \n'...
%                  'Using manually predefined set of channels\n'...
%                  'to reduce signal from neck area\n' ...
%                  '*********************************\n\n']);
%     fprintf(fidHTML,['<br><br>\n\n****  Detected use of HeadNeck_64 RF coil **** <br>\n'...
%                  'Using manually predefined set of channels<br>\n'...
%                  'to reduce signal from neck area<br>\n' ...
%                  '*********************************<br><br>\n\n']);    
% else
    iC_keep = 1:nc;
% end

nc_keep = length(iC_keep);
    



%% Check the number of FatNavs available compared to the size of the host data

nFatNavs = twix_obj.FatNav.dataSize(9)/twix_obj.image.NAve;

alignDim_rps = 2; % This is where it is assumed that the PE direction is the 'slow' acquisition direction
thisLine = squeeze(twix_obj.imageWithRefscan(1,1,:,1,1,1,1,1,1,1,1,1));
iSamp = find(thisLine);

fprintf(fidHTML,['No. of FatNavs: ' num2str(nFatNavs) '<br>\n']);
fprintf(fidHTML,['No. of measured lines in host sequence: ' num2str(length(iSamp)) '<br>\n']);
fprintf(['\n\nNo. of FatNavs: ' num2str(nFatNavs) '\n']);
fprintf(['No. of measured lines in host sequence: ' num2str(length(iSamp)) '\n']);

if length(iSamp)~=nFatNavs
    disp('Error: number of FatNavs doesn''t seem to match number of acquired lines found')
end

%% Process the FatNavs
% - First reconstruct each FatNav, then co-register using SPM to obtain
%   motion-estimates

[ACSims, timingReport_FatNavs, fatnavdir] = processFatNavs_GRAPPA4x4(twix_obj, ...
    outDir,'FatNavRes_mm',reconPars.FatNavRes_mm, 'iAve', reconPars.iAve, ...
    'appendString', appendString,'bSwapHandedness',reconPars.bSwapFatNavHandedness,...
    'bApplyNoseCircshift',reconPars.bApplyFatNavNoseCircshift,...
    'bUseNeckMasking',reconPars.bApplyFatNavNeckMasking,...
             'iAve', reconPars.iAve, 'appendString', appendString);

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

%% Do GRAPPA (if necessary) for the host sequence
% Note that this can use rather a lot of hard-disk space, as the raw data
% has to be shuffled around a few times to end up with one file per 
% reconstructed RF channel

if Arps(2) > 1
    useGRAPPAforHost = 1;
    if reconPars.bGRAPPAinRAM        
        [grappaRecon_1DFFT, mOutGRAPPA] = performHostGRAPPArecon2D_RAMonly(twix_obj,struct('iAve',reconPars.iAve,'iRep',reconPars.iRep));
        mOutGRAPPA.grappaRecon_1DFFT = grappaRecon_1DFFT; clear grappaRecon_1DFFT;
        timingReport_hostRecon = mOutGRAPPA.timingReport;
    else

        [timingReport_hostRecon, tempNameRoots] = performHostGRAPPArecon2D_toDisk(twix_obj,tempDir,struct('iAve',reconPars.iAve,'iRep',reconPars.iRep));    
    end
else
    useGRAPPAforHost = 0;
    timingReport_hostRecon = 0;
end



%% load moco parameters and align their orientation to the host data

fitResult = load([outDir '/motion_parameters_spm_' MIDstr appendString '.mat']);

this_fitMat = fitResult.MPos_cent.mats;

% Account for the fact that the host sequence may not have been acquired at
% isocentre:
%

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


%% Double-check the orientation with host data against the images from the FatNavs
% left/right matching
% - find the coil which has the biggest asymmetry left > right, and see if it
% is on the same side in both the FatNavs and the host data
nc_FatNavs = size(ACSims,4);
asymData_left = sum(reshape(abs(ACSims(1:floor(reconPars.FatNav_xyz(1)/2),:,:,:)),[],nc_FatNavs),1);
asymData_right = sum(reshape(abs(ACSims(ceil(reconPars.FatNav_xyz(1)/2):reconPars.FatNav_xyz(1),:,:,:)),[],nc_FatNavs),1);
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


%%
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



%% FatNavs are not acquired concurrently with the host data, so in the 
% % case of 'brittle' motion it may be especially beneficial to average
% % neighbouring motion estimates rather than taking the one from the same TR
 
rotTrans = rotmat2euler(this_fitMat_mm(1:3,1:3,:));
rotTrans(4:6,:) = squeeze(this_fitMat_mm(1:3,4,:));

newRotTrans = [rotTrans(:,1) (rotTrans(:,1:end-1)+rotTrans(:,2:end))/2];

this_fitMat_mm = zeros(size(this_fitMat_mm));
this_fitMat_mm(4,4,:) = 1;
    
this_fitMat_mm(1:3,1:3,:) = euler2rmat(newRotTrans(1:3,:));
this_fitMat_mm(1:3,4,:) = newRotTrans(4:6,:);
   


%% If GRAPPA was used in the 'slow' PE direction then the motion parameters will need to be interpolated to have values throughout k-space

%% THIS BIT COULD BE IMPROVED for handling 2D acceleration...

% 
% % if accelerated in fatnav direction, then motion needs to be interpolated
% % for the gaps
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



%% Apply the motion-correction by spawning jobs on the cluster
% (also need to spawn clean-up and finishing jobs which wait for this to
% finish)

RETROMOCOBOX_PATH = fileparts(which('addRetroMoCoBoxToPath.m'));
% mirtpath = fileparts(which('nufft.m'));
% MIRT_PATH = [mirtpath '/../'];
SPM_PATH = fileparts(which('spm.m')); 

tempNameRoots.allReconPars = [tempDir '/tempAllReconPars_' MIDstr '.mat'];
tempNameRoots.finalImage = [tempDir '/tempFinalImage_' MIDstr];
tempNameRoots.clusterScript = [tempDir '/tempClusterScript_' MIDstr '.sh'];
tempNameRoots.clusterScriptRecombine = [tempDir '/tempClusterScriptRecombine_' MIDstr '.sh'];
tempNameRoots.clusterScriptCleanup = [tempDir '/tempClusterScriptCleanup_' MIDstr '.sh'];
tempNameRoots.cleanupFiles = [tempDir '/tempCleanupPars_' MIDstr '.mat'];

iS = 1; % GRE data has only one 'set'
flipAxes_xyz = reconPars.flipAxes_xyz;
save(tempNameRoots.allReconPars,'iS','nc_keep','Hxyz','hxyz','nEco','reconPars','iC_keep','Hrps','hrps','permutedims',...
'fitMats_mm_toApply','alignDim_rps','alignIndices','hostVoxDim_mm','hostVoxDim_mm_rps', 'flipAxes_xyz','kspaceCentre_rps','tempNameRoots','MIDstr','outDir',...
'RETROMOCOBOX_PATH','SPM_PATH','tempDir')

fid = fopen(tempNameRoots.clusterScript,'w');
fprintf(fid,'#!/bin/bash\n');
fprintf(fid,['#SBATCH --array 1-' num2str(nEco) '\n']);
% fprintf(fid,['#SBATCH --array 1-' num2str(nEco) '%%1\n']); % use the %%1 to specify that only one can run at once for now to try to avoid parfor/job array problems
fprintf(fid,'#SBATCH -p cubric-centos7\n');
fprintf(fid,'#SBATCH --job-name=GREreconHelper\n');
fprintf(fid,['#SBATCH -o ' reconPars.CLUSTER_LOG_PATH '/GREreconHelperArray_%%A_%%a.out\n']);
fprintf(fid,['#SBATCH -e ' reconPars.CLUSTER_LOG_PATH '/GREreconHelperArray_%%A_%%a.err\n']);
fprintf(fid,'#SBATCH --ntasks 1\n'); 
% if nEco <= 10
%     fprintf(fid,['#SBATCH --cpus-per-task ' num2str(nEco) '\n']); % only 7 echoes, so parfor only needs 7 threads
% else
%     fprintf(fid,'#SBATCH --cpus-per-task 10\n'); 
% end
fprintf(fid,'#SBATCH --cpus-per-task 12\n'); 

%fprintf(fid,'#SBATCH --mem-per-cpu=32000M\n'); % MaxRSS showed requiring ~30 GB RAM for 336x336x192x7 data - I don't know why it still needs so much...
%%%% I even tried switching around the parfor loop in cluster_runMultieEchoGRE_eachcoil to try to get the RAM usage down - I'm not sure if this is a bug in this version of MATLAB...
fprintf(fid,'#SBATCH --mem=180G\n');  % 'seff' on jobid shows requiring around 27GB for 336x336x192x7 data or 80GB for 32 coils instead of 7 echoes
fprintf(fid,'#SBATCH --begin=now\n');
fprintf(fid,'#SBATCH --requeue\n');
fprintf(fid,['cd ' RETROMOCOBOX_PATH '/cluster\n']);
% fprintf(fid,['matlab -nodisplay -nodesktop -nosplash -r "cluster_runMultiEchoGRE_eachCoil(''' tempNameRoots.allReconPars ''',${SLURM_ARRAY_TASK_ID},1);exit;"\n']);
fprintf(fid,['matlab -nodisplay -nodesktop -nosplash -r "cluster_runMultiEchoGRE_eachEcho(''' tempNameRoots.allReconPars ''',${SLURM_ARRAY_TASK_ID});exit;"\n']);
fclose(fid);

disp('Launching batch job array for applying the motion correction...')
[status, sbatch_out] = system(['sbatch ' tempNameRoots.clusterScript]); 
jobnumber = sbatch_out(21:end); % after the text 'Submitted batch job ..'
disp('Done.')

fid = fopen(tempNameRoots.clusterScriptRecombine,'w');
fprintf(fid,'#!/bin/bash\n');
fprintf(fid,['#SBATCH --dependency=afterany:' jobnumber]); % wait for the array above to finish before starting the recombine
%fprintf(fid,'#SBATCH -p cubric-a100\n');
fprintf(fid,'#SBATCH -p cubric-rapids\n');
fprintf(fid,'#SBATCH --job-name=GREreconRecombine\n');
fprintf(fid,['#SBATCH -o ' reconPars.CLUSTER_LOG_PATH '/GREreconRecombine_%%j.out\n']);
fprintf(fid,['#SBATCH -e ' reconPars.CLUSTER_LOG_PATH '/GREreconRecombine_%%j.err\n']);
fprintf(fid,'#SBATCH --ntasks 1\n');
fprintf(fid,'#SBATCH --cpus-per-task 10\n');
%fprintf(fid,'#SBATCH --mem-per-cpu=32000M\n'); % this one also exceeded memory limit when set to 16000
%fprintf(fid,'#SBATCH --mem=320G\n');
fprintf(fid,'#SBATCH --mem=160G\n');    
% Here currently need at least enough to hold two full copies of full size
% recon - but 'seff' showed as much as 129GB for 336x336x192x7 data
fprintf(fid,'#SBATCH --begin=now\n');
fprintf(fid,['cd ' RETROMOCOBOX_PATH '/cluster\n']);
fprintf(fid,['matlab -nodisplay -nodesktop -nosplash -r "reconParsFile = ''' tempNameRoots.allReconPars ''';cluster_combineCoils_forASPIRE;exit;"\n']);
fclose(fid);

disp('Launching batch job array for applying the coil combination...')
[status, sbatch_out] = system(['sbatch ' tempNameRoots.clusterScriptRecombine]); 
jobnumber2 = sbatch_out(21:end); % after the text 'Submitted batch job ..'
disp('Done.')

ASPIRE_HOME = reconPars.ASPIRE_HOME;
save(tempNameRoots.cleanupFiles,'htmlDir','startTime','timingReport_FatNavs','nc_keep','nEco', ...
    'reconPars','tempNameRoots','tempDir','outDir','nFatNavs','fatnavdir','RETROMOCOBOX_PATH','ASPIRE_HOME','SPM_PATH');
fclose(fidHTML); % HTML needs closing and will be reopened in cleanup

% setenv('MRITOOLS_HOME','/home/scedg10/code/mritools_Linux_3.4.3/bin/')
fid = fopen(tempNameRoots.clusterScriptCleanup,'w');
fprintf(fid,'#!/bin/bash\n');
fprintf(fid,['#SBATCH --dependency=afterany:' jobnumber2]); % wait for the recombine to finish before starting the cleanup
fprintf(fid,'#SBATCH -p cubric-centos7\n');
fprintf(fid,'#SBATCH --job-name=GREreconCleanup\n');
fprintf(fid,['#SBATCH -o ' reconPars.CLUSTER_LOG_PATH '/GREreconCleanup_%%j.out\n']);
fprintf(fid,['#SBATCH -e ' reconPars.CLUSTER_LOG_PATH '/GREreconCleanup_%%j.err\n']);
fprintf(fid,'#SBATCH --ntasks 1\n');
fprintf(fid,'#SBATCH --mem-per-cpu=32G\n');
fprintf(fid,'#SBATCH --begin=now\n');
fprintf(fid,['cd ' RETROMOCOBOX_PATH '/cluster\n']);
fprintf(fid,['matlab -nodisplay -nodesktop -nosplash -r "cleanupFile = ''' tempNameRoots.cleanupFiles ''';cluster_cleanup_andASPIRE;exit;"\n']);
fclose(fid);

disp('Launching batch job array for applying the cleanup...')
[status, sbatch_out] = system(['sbatch ' tempNameRoots.clusterScriptCleanup]); 
disp('Done.')
