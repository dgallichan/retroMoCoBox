function timingReport = reconstructSiemensVolume(twix_obj,reconPars)
% function timingReport = reconstructSiemensVolume(twix_obj,reconPars)
% 
% Called by reconstructSiemensMP2RAGEwithFatNavs.m to reconstruct a
% specific 'average' or 'repetition'.
%
% Feb 2018, gallichand@cardiff.ac.uk



retroMocoBoxVersion = reconPars.retroMocoBoxVersion; % put this into the HTML for reference

%%

if ~isfield(reconPars,'iAve')
    reconPars.iAve = 1;
end
if ~isfield(reconPars,'iRep')
    reconPars.iRep = 1;
end


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

figIndex = 999; % figure to use for outputs

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

outDir = [reconPars.outRoot '/hostrecon_' MIDstr appendString];
if ~exist(outDir,'dir')
    mkdir(outDir)
end

htmlDir = [outDir '/html'];
if ~exist(htmlDir,'dir')
    mkdir(htmlDir)
end

if isempty(reconPars.tempRoot)
    reconPars.tempRoot = reconPars.outRoot;
end
tempDir = [reconPars.tempRoot '/temp_' MIDstr appendString];
if ~exist(tempDir,'dir')
    mkdir(tempDir)
end

% intialize html index file
fid = fopen([htmlDir '/index.html'],'w');
fprintf(fid,['<html><head><title>MP(2)RAGE with FatNavs - Summary</title>\n']);
fprintf(fid,'</head>\n');
fprintf(fid,'<body>\n');
fprintf(fid,['<h2>MP(2)RAGE with FatNavs - Summary: %s</h2>\n'],[fileName]);
fprintf(fid,['<br>\n']);





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
        FatNav_xyz = FatNav_FOVxyz ./ reconPars.FatNavRes_mm;
    case 6
        FatNav_FOVxyz = [192 264 264]; % FatNav FOV         
        FatNav_xyz = FatNav_FOVxyz ./ reconPars.FatNavRes_mm;
        FatNav_xyz(3) = 64; % No idea why, but the 6mm data has 64 points in the readout direction for the ACS lines instead of 44...        
end



%%


% in the next line I assume that Siemens always use the same structure for
% the 'tReferenceImage0' field - but I haven't looked for documentation to
% support this, so it may not always extract the scan date properly...
if isfield(twix_obj.hdr.MeasYaps,'tReferenceImage0') % some files apparently might not even have this field...
    iDot = strfind(char(twix_obj.hdr.MeasYaps.tReferenceImage0),'.'); % some versions seem to create a char already, others not...
    thisDateString = char(twix_obj.hdr.MeasYaps.tReferenceImage0);
    if strcmp(thisDateString(iDot(end)+(1:6)),'300000') % it seems sometimes this 300000 is present, then the year is only two digits...
        fprintf(fid,['<strong>Date of scan:</strong> ' local_reformatDateString(['20' thisDateString(iDot(end)+6+(1:6))]) '<br>\n']);
    else
        fprintf(fid,['<strong>Date of scan:</strong> ' local_reformatDateString(thisDateString(iDot(end)+(1:8))) '<br>\n']);
    end
end

if ischar(reconPars.bKeepPatientInfo)
    fprintf(fid,['<strong>Database ID:</strong> ' reconPars.bKeepPatientInfo '<br>\n']);
else
    if reconPars.bKeepPatientInfo
        %% Put patient info into HTML
        
      if isfield(twix_obj.hdr.Config,'PatientName')
            fprintf(fid,['<strong>Patient Name:</strong> ' twix_obj.hdr.Config.PatientName '<br>\n']);
        else
            if isfield(twix_obj.hdr.Config,'tPatientName') % why would this be different for some scanners / software versions...?!?
                fprintf(fid,['<strong>Patient Name:</strong> ' twix_obj.hdr.Config.tPatientName '<br>\n']);
            end
        end
        fprintf(fid,['<strong>Patient ID:</strong> ' twix_obj.hdr.Config.PatientID '<br>\n']);
        switch twix_obj.hdr.Config.PatientSex
            case 1
                fprintf(fid,['<strong>Patient Sex:</strong> Female<br>\n']);
            case 2
                fprintf(fid,['<strong>Patient Sex:</strong> Male<br>\n']);
            case 3
                fprintf(fid,['<strong>Patient Sex:</strong> Other<br>\n']);
        end
        fprintf(fid,['<strong>Patient date of birth:</strong> ' local_reformatDateString(num2str(twix_obj.hdr.Config.PatientBirthDay)) '<br>\n']);
        
    else
        fprintf(fid,['<strong>Patient Info:</strong> Anonymised<br>\n']);
    end
end


%% Try to establish slice orientation - this has not been rigourously tested for all possible orientations...!

sNormal = twix_obj.hdr.MeasYaps.sSliceArray.asSlice{1}.sNormal;
if ~isfield(sNormal,'dCor'),  sNormal.dCor = 0; end
if ~isfield(sNormal,'dTra'),  sNormal.dTra = 0; end
if ~isfield(sNormal,'dSag'),  sNormal.dSag = 0; end

allNormals = [sNormal.dCor, sNormal.dTra, sNormal.dSag].';

baseOrientation = find(allNormals==max(allNormals));% 1 - Coronal; 2 - Transverse, 3 = Sagittal

if ~isfield(twix_obj.hdr.MeasYaps.sSliceArray.asSlice{1},'sPosition')
    sPosition.dCor = 0;
    sPosition.dTra = 0;
    sPosition.dSag = 0;
else
    sPosition = twix_obj.hdr.MeasYaps.sSliceArray.asSlice{1}.sPosition;
    if ~isfield(sPosition,'dCor'), sPosition.dCor = 0; end
    if ~isfield(sPosition,'dTra'), sPosition.dTra = 0; end
    if ~isfield(sPosition,'dSag'), sPosition.dSag = 0; end
end


permutedims = []; % initialize
if isfield(twix_obj.hdr.MeasYaps.sSliceArray.asSlice{1},'dInPlaneRot')
    dInPlaneRot = twix_obj.hdr.MeasYaps.sSliceArray.asSlice{1}.dInPlaneRot;
else
    dInPlaneRot = 0;
end
switch baseOrientation 
    case 1 % Coronal
        orientText = 'coronal';
        if (dInPlaneRot*180/pi) < 45
            permutedims = [2 3 1];
        else
            permutedims = [1 3 2]; 
        end       
    case 2 % Transverse
        orientText = 'transversal';
        if (dInPlaneRot*180/pi) < 45
            permutedims = [1 2 3];
        else
            permutedims = [2 1 3];
        end
    case 3 % Sagittal       
        orientText = 'sagittal';
        if (dInPlaneRot*180/pi) < 45
            permutedims = [3 2 1];
        else
            permutedims = [3 1 2];
        end
end

% reconPars.swapDims_xyz = [0 0 1]; % <-- so far all data I have seems to fit this - only the z direction needs reversal - 
                        %     which would presumably be equivalent to the Siemens definition of z as running from head to foot
                        % 4/1/16 - now specify this as an input because I
                        % have found some data which didn't fit this
                        % default...
[~, permutedims_toXYZ] = sort(permutedims);
swapDims_rps = reconPars.swapDims_xyz(permutedims_toXYZ);

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

kspaceCentre_rps = [ceil(twix_obj.image.centerCol(1)/2) twix_obj.image.centerLin(1) twix_obj.image.centerPar(1)];
for iD = 1:3
    if swapDims_rps(iD)
        kspaceCentre_rps(iD) = hrps(iD)-kspaceCentre_rps(iD);
    end
end

kspaceCentre_xyz = kspaceCentre_rps(permutedims);

nc = twix_obj.image.NCha;
nS = twix_obj.image.NSet; % for MP2RAGE this will be 2...


% and put stuff into html report and the screen
if nS == 2
    fprintf(fid,['<h4>Host MP2RAGE ']);
    fprintf(['\n\n\nHost MP2RAGE ']);
else
    fprintf(fid,['<h4>Host MPRAGE ']);
    fprintf(['\n\n\nHost MPRAGE ']);
end
fprintf(fid,['sequence</h4>\n']);
fprintf(['sequence\n']);

fprintf(fid,['RPS dimensions in raw data: ' num2str(hrps(1)) 'x' num2str(hrps(2)) 'x' num2str(hrps(3)) '<br>\n']);
fprintf(['RPS dimensions in raw data: ' num2str(hrps(1)) 'x' num2str(hrps(2)) 'x' num2str(hrps(3)) '\n']);
fprintf(fid,['XYZ dimensions reconstructed: ' num2str(Hxyz(1)) 'x' num2str(Hxyz(2)) 'x' num2str(Hxyz(3)) '<br>\n']);
fprintf(['XYZ dimensions reconstructed: ' num2str(Hxyz(1)) 'x' num2str(Hxyz(2)) 'x' num2str(Hxyz(3)) '\n']);
fprintf(fid,['FOV - ' num2str(FOVxyz(1),'%.1f') 'x' num2str(FOVxyz(2),'%.1f') 'x' num2str(FOVxyz(3),'%.1f') 'mm<br>\n']);
fprintf(['FOV - ' num2str(FOVxyz(1),'%.1f') 'x' num2str(FOVxyz(2),'%.1f') 'x' num2str(FOVxyz(3),'%.1f') 'mm\n']);
fprintf(fid,['Resolution: ' num2str(hostVoxDim_mm(1),'%.3f') 'x' num2str(hostVoxDim_mm(2),'%.3f') 'x' num2str(hostVoxDim_mm(3),'%.3f') 'mm<br>\n']);
fprintf(['Resolution: ' num2str(hostVoxDim_mm(1),'%.3f') 'x' num2str(hostVoxDim_mm(2),'%.3f') 'x' num2str(hostVoxDim_mm(3),'%.3f') 'mm\n']);
if manualFatNavRes
    fprintf(fid,['Manually selected FatNav resolution: ' num2str(reconPars.FatNavRes_mm) ' mm<br>\n']);
    fprintf(['Manually selected FatNav resolution: ' num2str(reconPars.FatNavRes_mm) ' mm\n']);
else
    fprintf(fid,['Assumed FatNav resolution (based on field strength): ' num2str(reconPars.FatNavRes_mm) ' mm<br>\n']);
    fprintf(['Assumed FatNav resolution (based on field strength): ' num2str(reconPars.FatNavRes_mm) ' mm\n']);
end
fprintf(fid,['Detected orientation: ' orientText '<br>\n']);
fprintf(['Detected orientation: ' orientText '\n']);
fprintf(fid,['SwapDims_xyz: [' num2str(reconPars.swapDims_xyz(1)) ' ' num2str(reconPars.swapDims_xyz(2)) ' ' num2str(reconPars.swapDims_xyz(3)) ']<br>\n']); 
fprintf(['SwapDims_xyz: [' num2str(reconPars.swapDims_xyz(1)) ' ' num2str(reconPars.swapDims_xyz(2)) ' ' num2str(reconPars.swapDims_xyz(3)) ']\n']); 
if isfield(twix_obj.hdr.MeasYaps,'sCoilSelectMeas')
    fprintf(['Coil used: ' char(twix_obj.hdr.MeasYaps.sCoilSelectMeas.aRxCoilSelectData{1}.asList{1}.sCoilElementID.tCoilID) ', with ' num2str(nc) ' channels active\n']);
    fprintf(fid,['Coil used: ' char(twix_obj.hdr.MeasYaps.sCoilSelectMeas.aRxCoilSelectData{1}.asList{1}.sCoilElementID.tCoilID) ', with ' num2str(nc) ' channels active<br>\n']);
else
    fprintf(['Coil used: ' char(twix_obj.hdr.MeasYaps.asCoilSelectMeas{1}.asList{1}.sCoilElementID.tCoilID) ', with ' num2str(nc) ' channels active\n']);
    fprintf(fid,['Coil used: ' char(twix_obj.hdr.MeasYaps.asCoilSelectMeas{1}.asList{1}.sCoilElementID.tCoilID) ', with ' num2str(nc) ' channels active<br>\n']);
end

%% Check if using HEADNECK_64 receive coil, and discard channels over the neck if this is the case (would be nice to know what Siemens does...)
if isfield(twix_obj.hdr.MeasYaps,'sCoilSelectMeas') ...
        && strcmp(twix_obj.hdr.MeasYaps.sCoilSelectMeas.aRxCoilSelectData{1}.asList{1}.sCoilElementID.tCoilID,'"HeadNeck_64"')
    
    iC_keep = 1:nc;
    if nc==52
        iC_keep([1 7 8 18 29 30 39 40 49 50]) = []; % these channels cover the neck (in one test dataset with 52 data channels from the 64-channel coil...)
        fprintf(['\n\n****  Detected use of HeadNeck_64 RF coil **** \n'...
            'Using manually predefined set of channels\n'...
            'to reduce signal from neck area\n' ...
            '*********************************\n\n']);
        fprintf(fid,['<br><br>\n\n****  Detected use of HeadNeck_64 RF coil **** <br>\n'...
            'Using manually predefined set of channels<br>\n'...
            'to reduce signal from neck area<br>\n' ...
            '*********************************<br><br>\n\n']);
    end
else
    iC_keep = 1:nc;
end

nc_keep = length(iC_keep);
    




%% Check the number of FatNavs available compared to the size of the host data

% find the right (chronological) order for the k-space lines in the raw data

if reconPars.bLinParSwap % MP2RAGE sequence currently doesn't allow 2D GRAPPA, so this means the number of FatNavs should already match that dimension
    alignDim = permutedims_toXYZ(3); 
    iSamp = 1:hxyz(alignDim);           
else
    alignDim = permutedims_toXYZ(2);    
%     iSamp = twix_obj.imageWithRefscan.Lin(1:hrps(3)*nS:end);
    thisLine = squeeze(twix_obj.imageWithRefscan(1,1,:,1,1,1,1,1,1,1,1,1));
    iSamp = find(thisLine);
end
nFatNavs = twix_obj.FatNav.dataSize(9)/twix_obj.image.NAve;

fprintf(fid,['No. of FatNavs: ' num2str(nFatNavs) '<br>\n']);
fprintf(fid,['No. of measured lines in host sequence: ' num2str(length(iSamp)) '<br>\n']);
fprintf(['\n\nNo. of FatNavs: ' num2str(nFatNavs) '\n']);
fprintf(['No. of measured lines in host sequence: ' num2str(length(iSamp)) '\n']);

if length(iSamp)~=nFatNavs
    disp('Error: number of FatNavs doesn''t seem to match number of acquired lines found')
end

%% Process the FatNavs
% - First reconstruct each FatNav, then co-register using SPM to obtain
%   motion-estimates

[ACSims, timingReport_FatNavs, fatnavdir] = processFatNavs_GRAPPA4x4(twix_obj, ... 
            reconPars.outRoot,'FatNavRes_mm',reconPars.FatNavRes_mm, 'iAve', reconPars.iAve, 'appendString', appendString);

% And put stuff into the html report
if exist(fatnavdir,'dir') % could have just kept the motion-parameters file...
    imdims = 2*[304 128];
    
    if exist([fatnavdir '/a_FatNav_ACSim_' MIDstr '.png'],'file')
        copyfile([fatnavdir '/a_FatNav_ACSim_' MIDstr '.png'],[htmlDir '/ACSim.png']);
        fprintf(fid,['ACS image:<br>\n']);
        fprintf(fid,['<img src="ACSim.png" height=%s width=%s><br><br>\n'],num2str(imdims(2)),num2str(imdims(1)));
    end
    if exist([fatnavdir '/a_FatNav_MoCoPars_' MIDstr '.png'],'file')
        copyfile([fatnavdir '/a_FatNav_MoCoPars_' MIDstr '.png'],[htmlDir '/motion_parameters.png']);
        fprintf(fid,['<h4>FatNavs</h4>\n']);
        fprintf(fid,['Estimated motion parameters:<br>\n']);
        fprintf(fid,['<img src="motion_parameters.png"><br><br>\n']);
    end
    if exist([fatnavdir '/a_mov_eachFatNav.gif'],'file')
        copyfile([fatnavdir '/a_mov_eachFatNav.gif'],[htmlDir '/mov_eachFatNav.gif'])
        fprintf(fid,['15 example FatNavs covering complete scan:<br>\n']);
        fprintf(fid,['<img src="mov_eachFatNav.gif" height=%s width=%s><br><br>\n'],num2str(imdims(2)),num2str(imdims(1)));
    end
    if exist([fatnavdir '/a_mov_spm_eachFatNav.gif'],'file')
        copyfile([fatnavdir '/a_mov_spm_eachFatNav.gif'],[htmlDir '/mov_spm_eachFatNav.gif'])
        fprintf(fid,['And after rigid-body registration using SPM:<br>\n']);
        fprintf(fid,['<img src="mov_spm_eachFatNav.gif" height=%s width=%s><br><br>\n'],num2str(imdims(2)),num2str(imdims(1)));
    end
end

%% Do GRAPPA (if necessary) for the host sequence
% Note that this can use rather a lot of hard-disk space, as the raw data
% has to be shuffled around a few times to end up with one file per 
% reconstructed RF channel

if Arps(2) > 1
    useGRAPPAforHost = 1;
    if reconPars.bGRAPPAinRAM        
        [grappaRecon_1DFFT, mOutGRAPPA] = performHostGRAPPArecon_RAMonly(twix_obj,struct('iAve',reconPars.iAve,'iRep',reconPars.iRep));
        mOutGRAPPA.grappaRecon_1DFFT = grappaRecon_1DFFT; clear grappaRecon_1DFFT;
        timingReport_hostRecon = mOutGRAPPA.timingReport;
    else
%         [mOutGRAPPA, timingReport_hostRecon] =
%         performHostGRAPPArecon(twix_obj,tempDir); % the original low-RAM
%         alternative code needs verifying as the output ends up corrupted
%         (in at least some cases)
        [timingReport_hostRecon, tempNameRoots] = performHostGRAPPArecon_toDisk(twix_obj,tempDir,struct('iAve',reconPars.iAve,'iRep',reconPars.iRep));    
    end
else
    useGRAPPAforHost = 0;
    timingReport_hostRecon = 0;
end


%% load moco parameters and align their orientation to the host data

fitResult = load([reconPars.outRoot '/motion_parameters_spm_' MIDstr appendString '.mat']);

this_fitMat = fitResult.MPos_cent.mats;

% Account for the fact that the host sequence may not have been acquired at
% isocentre:
%
%%% NOTE - I have not checked the following code works properly in the case
%%% of non-isotropic host FOV and/or voxels situated away from isocentre.
vox_fitMat = this_fitMat;
vox_fitMat(1:3,4,:) = bsxfun(@times,this_fitMat(1:3,4,:), 1./hostVoxDim_mm); % rescale translations into new voxel size

rotMatDisplacement_mm = [sPosition.dSag sPosition.dCor -sPosition.dTra].';
rotMatDisplacement_voxels = rotMatDisplacement_mm ./ hostVoxDim_mm;

new_fitMat = zeros(size(vox_fitMat));
for iT = 1:size(vox_fitMat,3)
    thisRotMat = vox_fitMat(1:3,1:3,iT);
    thisRotMatDisplacement = rotMatDisplacement_voxels - (thisRotMat*rotMatDisplacement_voxels);
    new_fitMat(1:3,1:3,iT) = thisRotMat;
    new_fitMat(1:3,4,iT) = thisRotMatDisplacement + vox_fitMat(1:3,4,iT);
end

this_fitMat_mm = new_fitMat;
this_fitMat_mm(1:3,4,:) = bsxfun(@times,this_fitMat_mm(1:3,4,:), hostVoxDim_mm); % put translations into mm

theseDisplacements = squeeze(this_fitMat_mm(1:3,4,:));

switch baseOrientation
    case 1 % Coronal         
        tiltTheta = acos(allNormals(1))*180/pi;
        thisRot = euler2rmat(tiltTheta,0,0);        
    case 2 % Transverse
        tiltTheta = acos(allNormals(2))*180/pi;
        thisRot = euler2rmat(tiltTheta,0,0);
    case 3 % Sagittal
        thisRot = euler2rmat(dInPlaneRot*180/pi,0,0);
end

% and account for rotation of slices...
newDisplacements = thisRot*theseDisplacements;
this_fitMat_mm(1:3,4,:) = newDisplacements;



%% Double-check the orientation with host data against the images from the FatNavs
% left/right matching
% - find the coil which has the biggest asymmetry left > right, and see if it
% is on the same side in both the FatNavs and the host data
nc_FatNavs = size(ACSims,4);
asymData_left = sum(reshape(abs(ACSims(1:floor(FatNav_xyz(1)/2),:,:,:)),[],nc_FatNavs),1);
asymData_right = sum(reshape(abs(ACSims(ceil(FatNav_xyz(1)/2):FatNav_xyz(1),:,:,:)),[],nc_FatNavs),1);
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
hostExampleVolume = permute(hostExampleVolume,permutedims);

if reconPars.swapDims_xyz(1)
    hostExampleVolume = hostExampleVolume(end:-1:1,:,:);
end

%%
ov1 = orthoview(ACSims(:,:,:,iAsymCoil),'drawIms',0,'mip',1,'clims',[0 max(abs(ACSims(:)))]);
ov2 = orthoview(hostExampleVolume,'drawIms',0,'mip',0,'clims',[0 max(abs(hostExampleVolume(:)))]);
fig(figIndex)
clf
set(gcf,'Position',[    22   594   702   473])
subplot1(1,2)
subplot1(1)
imab(ov1.im1)
title(['xy MIP of FatNav, coil ' num2str(iAsymCoil)])
subplot1(2)
imab(ov2.im1)
title(['xy centre-slice of host GRAPPA recon, coil ' num2str(iAsymCoil)])
colormap(gray)
export_fig([htmlDir '/orientationCheck_xy.png']);
fprintf(fid,['Orientation check for left/right symmetry:<br>\n']);
fprintf(fid,['<img src="orientationCheck_xy.png"><br><br>\n']);
fprintf(fid,['(Both images above should have the brightest signal on the left of the image. If not, the orientation of the FatNavs is not correctly aligned with the host sequence)<br><br><br>\n']);


%%

if exist([fatnavdir '/eachFatNav_001.nii'],'file')
    affMat = thisRot;
%     affMat(1:3,4) = rotMatDisplacement_mm.'/2; % Why was this divided by
%     two...?
    affMat(1:3,4) = rotMatDisplacement_mm.'; 
    affMat(4,4) = 1;
    
    fileTest1 = [fatnavdir '/test1.nii'];
    fileTest2 = [fatnavdir '/test2.nii'];
    copyfile([fatnavdir '/eachFatNav_001.nii'],fileTest1);
    copyfile([fatnavdir '/eachFatNav_001.nii'],fileTest2);
    
    V = spm_vol_nifti(fileTest1);
    V(2) = spm_vol_nifti(fileTest2);
    newMat = V(1).mat;
    newMat = affMat*newMat;
    V(2).mat = newMat;
    spm_reslice(V,struct('mask',false,'mean',false,'interp',5,'which',1,'wrap',[0 0 0],'prefix',''));
    newFatIm = rn(fileTest2);
    
    nShowX = min(FatNav_xyz(1),round(0.5*FatNav_xyz(1)*FOVxyz(1)/FatNav_FOVxyz(1))*2);
    nShowY = min(FatNav_xyz(2),round(0.5*FatNav_xyz(2)*FOVxyz(2)/FatNav_FOVxyz(2))*2);
    nShowZ = min(FatNav_xyz(3),round(0.5*FatNav_xyz(3)*FOVxyz(3)/FatNav_FOVxyz(3))*2);
    
    fig(figIndex)
    orthoview(newFatIm([1:nShowX]-nShowX/2+FatNav_xyz(1)/2,[1:nShowY]-nShowY/2+FatNav_xyz(2)/2,[1:nShowZ]-nShowZ/2+FatNav_xyz(3)/2),'useNewFig',0);
    set(gcf,'Position',[    50   720   950  340])
    subplot1(2)
    title({'First FatNav realigned to coordinate system of host sequence','xz'})
    export_fig([htmlDir '/orientationCheck_FatVolume.png'])
    fprintf(fid,['Orientation check for host sequence slice rotation and positioning:<br>\n']);
    fprintf(fid,['<img src="orientationCheck_FatVolume.png"><br><br>\n']);
    fprintf(fid,['(The fat volume shown above should approximately correspond to the FOV chosen for the host sequence)<br><br><br>\n']);
end

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
if nFatNavs ~= hxyz(alignDim)
    rotTrans = rotmat2euler(this_fitMat_mm(1:3,1:3,:));
    rotTrans(4:6,:) = squeeze(this_fitMat_mm(1:3,4,:));
    
    rotTrans_interp = interp1(iSamp,rotTrans.',1:hxyz(alignDim)).';
    startIndex = iSamp(1);
    % copy in values at edges of k-space to avoid extrapolation errors
    if startIndex > 1
        rotTrans_interp(:,1:startIndex-1) = repmat(rotTrans_interp(:,startIndex),[1 startIndex-1]);
    end
    endIndex = iSamp(end);
    if endIndex < hxyz(alignDim)
        rotTrans_interp(:,endIndex+1:end) = repmat(rotTrans_interp(:,endIndex),[1 hxyz(alignDim)-endIndex]);
    end
    
    this_fitMat_mm_interp = zeros(4,4,hxyz(alignDim));
    this_fitMat_mm_interp(4,4,:) = 1;
    
    this_fitMat_mm_interp(1:3,1:3,:) = euler2rmat(rotTrans_interp(1:3,:));
    this_fitMat_mm_interp(1:3,4,:) = rotTrans_interp(4:6,:);
    
    fitMats_mm_toApply = this_fitMat_mm_interp;
    
else
    fitMats_mm_toApply = this_fitMat_mm;
    fitMats_mm_toApply(4,4,:) = 1;
end


if reconPars.swapDims_xyz(alignDim)
    % make it so that centre of k-space is not 'moved' (accounting for partial Fourier):
    fitMats_mm_toApply = recentre_affmats(fitMats_mm_toApply,hxyz(alignDim)-kspaceCentre_xyz(alignDim));  
    alignIndices = hxyz(alignDim):-1:1;  
else
    % make it so that centre of k-space is not 'moved' (accounting for partial Fourier):
    fitMats_mm_toApply = recentre_affmats(fitMats_mm_toApply,kspaceCentre_xyz(alignDim));  
    alignIndices = 1:hxyz(alignDim);
end



%% Prepare for the retrospective motion-correction of the host sequence
tStart_applyMoco = clock;

timingReport_applyMoco = zeros(nc,nS); % store the time to reconstruct each volume with the NUFFT

fprintf('............\n')
fprintf('... Performing NUFFT on each volume to apply motion correction parameters\n');


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

    for iC = 1:nc_keep           
        
        mOut.thiscoil_ims = complex(zeros(Hxyz(1),Hxyz(2),Hxyz(3),nS,'single'));
        mOut.thiscoil_ims_corrected = complex(zeros(Hxyz(1),Hxyz(2),Hxyz(3),nS,'single'));
        
        for iS = 1:nS
        
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
                % thisData = squeeze(mOutGRAPPA.dataCombined(:,:,:,2)); iC=1;iS = 1; % use this to have more brain coverage for debugging...
                thisData = fft1s(thisData,1); % put into full 3D k-space
            else
                thisData = squeeze(twix_obj.image(:,iC_keep(iC),:,:,1,reconPars.iAve,1,1,reconPars.iRep,iS));
            end
            
            thisData = permute(thisData,permutedims);
            if reconPars.swapDims_xyz(1)
                thisData = thisData(end:-1:1,:,:,:);
            end
            if reconPars.swapDims_xyz(2)
                thisData = thisData(:,end:-1:1,:,:);
            end
            if reconPars.swapDims_xyz(3)
                thisData = thisData(:,:,end:-1:1,:);
            end
            
            newData = thisData;
            
            if any(hxyz~=Hxyz)
                thisData(Hxyz(1),Hxyz(2),Hxyz(3)) = 0; % extend to new size
                thisData = circshift(thisData,double([Hxyz(1)-hxyz(1) Hxyz(2)-hxyz(2) Hxyz(3)-hxyz(3)].*(1-reconPars.swapDims_xyz)));
                % the 'double' in the above line appears necessary in certain
                % versions of Matlab, no idea why...
            end
            
            thisData = ifft3s(thisData)*prod(hxyz);
            
            tic
            thisData_corrected = applyRetroMC_nufft(newData,fitMats_mm_toApply,alignDim,alignIndices,11,hostVoxDim_mm,Hxyz,kspaceCentre_xyz);
            % the multi-CPU (parfor) version of the NUFFT can be called with a
            % negative value for the 'useTable' input:
            %    thisData_corrected = applyRetroMC_nufft(newData,fitMats_mm_toApply,alignDim,alignIndices,-11,hostVoxDim_mm,Hxyz,kspaceCentre_xyz);
            % However, in my latest tests (September 2016) I found that it was
            % actually slower for the current NUFFT parameters than the
            % non-parallelized version, so I've disabled it again as a
            % default...
            timingReport_applyMoco(iC,iS) = toc;
            
            if nS > 1
                mOut.thiscoil_ims(:,:,:,iS) = thisData;
                mOut.thiscoil_ims_corrected(:,:,:,iS) = thisData_corrected;
            else
                mOut.thiscoil_ims = thisData;
                mOut.thiscoil_ims_corrected = thisData_corrected;
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
        mOut.all_uniImage_corrected = mOut.all_uniImage_corrected./mOut.all_refImage_corrected;
        
        sn(int16(4095*mOut.all_ims(:,:,:,1)/max(reshape(mOut.all_ims,[],1))),[outDir '/INV1'],hostVoxDim_mm)
        sn(int16(4095*mOut.all_ims_corrected(:,:,:,1)/max(reshape(mOut.all_ims_corrected,[],1))),[outDir '/INV1_corrected'],hostVoxDim_mm)
        
        
        sn(int16(4095*(mOut.all_uniImage+0.5)) ,[outDir '/UNI'],hostVoxDim_mm)
        sn(int16(4095*(mOut.all_uniImage_corrected+0.5)),[outDir '/UNI_corrected'],hostVoxDim_mm)
        sn(int16(4095*mOut.all_ims(:,:,:,2)/max(reshape(mOut.all_ims,[],1))),[outDir '/INV2'],hostVoxDim_mm)
        sn(int16(4095*mOut.all_ims_corrected(:,:,:,2)/max(reshape(mOut.all_ims_corrected,[],1))),[outDir '/INV2_corrected'],hostVoxDim_mm)
    else
        % if using matfile variables you can't use index in dimensions which
        % aren't there, even if you only put a 1 there...!
        sn(int16(4095*mOut.all_ims/max(reshape(mOut.all_ims,[],1))),[outDir '/INV1'],hostVoxDim_mm)
        sn(int16(4095*mOut.all_ims_corrected/max(reshape(mOut.all_ims_corrected,[],1))),[outDir '/INV1_corrected'],hostVoxDim_mm)
        
    end
    
    
else % the much faster version with much hungrier RAM requirements:
    
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
            
            thisData = permute(thisData,permutedims);
            if reconPars.swapDims_xyz(1)
                thisData = thisData(end:-1:1,:,:,:);
            end
            if reconPars.swapDims_xyz(2)
                thisData = thisData(:,end:-1:1,:,:);
            end
            if reconPars.swapDims_xyz(3)
                thisData = thisData(:,:,end:-1:1,:);
            end
            
            newData = thisData;
            
            if any(hxyz~=Hxyz)
                thisData(Hxyz(1),Hxyz(2),Hxyz(3)) = 0; % extend to new size
                thisData = circshift(thisData,double([Hxyz(1)-hxyz(1) Hxyz(2)-hxyz(2) Hxyz(3)-hxyz(3)].*(1-reconPars.swapDims_xyz)));
                % the 'double' in the above line appears necessary in certain
                % versions of Matlab, no idea why...
            end
            
            thisData = ifft3s(thisData)*prod(hxyz);
            
            tic
            thisData_corrected = applyRetroMC_nufft(newData,fitMats_mm_toApply,alignDim,alignIndices,11,hostVoxDim_mm,Hxyz,kspaceCentre_xyz);
            timingReport_applyMoco(iC,iS) = toc;
            
            if nS > 1
                thiscoil_ims(:,:,:,iS) = thisData;
                thiscoil_ims_corrected(:,:,:,iS) = thisData_corrected;
            else
                thiscoil_ims = thisData;
                thiscoil_ims_corrected = thisData_corrected;
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
    all_uniImage_corrected = all_uniImage_corrected./all_refImage_corrected;
    
    sn(int16(4095*all_ims(:,:,:,1)/max(reshape(all_ims,[],1))),[outDir '/INV1'],hostVoxDim_mm)
    sn(int16(4095*all_ims_corrected(:,:,:,1)/max(reshape(all_ims_corrected,[],1))),[outDir '/INV1_corrected'],hostVoxDim_mm)
    
    if nS>1
        sn( int16(4095*(all_uniImage+0.5)) ,[outDir '/UNI'],hostVoxDim_mm)
        sn( int16(4095*(all_uniImage_corrected+0.5)),[outDir '/UNI_corrected'],hostVoxDim_mm)
        sn(int16(4095*all_ims(:,:,:,2)/max(reshape(all_ims,[],1))),[outDir '/INV2'],hostVoxDim_mm)
        sn(int16(4095*all_ims_corrected(:,:,:,2)/max(reshape(all_ims_corrected,[],1))),[outDir '/INV2_corrected'],hostVoxDim_mm)
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

switch nS
    case 1 % again different code because matfile variables can't index dimensions which don't exist
        
        ov1 = orthoview(mOut.all_ims,'drawIms',0);
        imab_overwrite([htmlDir '/INV1.png'],ov1.oneIm);
        ov1 = orthoview(mOut.all_ims_corrected,'drawIms',0);
        imab_overwrite([htmlDir '/INV1_corrected.png'],ov1.oneIm);
        
        fprintf(fid,['INV1 image before correction:<br>\n']);
        fprintf(fid,['<img src="INV1.png"><br><br>\n']);
        fprintf(fid,['INV1 image after correction:<br>\n']);
        fprintf(fid,['<img src="INV1_corrected.png"><br><br>\n']);
        
        
        testMagick = system('convert -version');
        
        if testMagick==0 % can use ImageMagick to make animated GIFs...
            processString = ['convert -dispose 2 -delay 50 -loop 0 ' htmlDir '/INV1*.png ' htmlDir '/mov_INV1.gif'];
            system(processString);
            fprintf(fid,['INV1 image movie before/after correction:<br>\n']);
            fprintf(fid,['<img src="mov_INV1.gif"><br><br>\n']);
        end
        
    case 2
        
        ov1 = orthoview(mOut.all_ims(:,:,:,1),'drawIms',0);
        imab_overwrite([htmlDir '/INV1.png'],ov1.oneIm);
        ov1 = orthoview(mOut.all_ims_corrected(:,:,:,1),'drawIms',0);
        imab_overwrite([htmlDir '/INV1_corrected.png'],ov1.oneIm);
        
        fprintf(fid,['INV1 image before correction:<br>\n']);
        fprintf(fid,['<img src="INV1.png"><br><br>\n']);
        fprintf(fid,['INV1 image after correction:<br>\n']);
        fprintf(fid,['<img src="INV1_corrected.png"><br><br>\n']);
        
        
        testMagick = system('convert -version');
        
        if testMagick==0 % can use ImageMagick to make animated GIFs...
            processString = ['convert -dispose 2 -delay 50 -loop 0 ' htmlDir '/INV1*.png ' htmlDir '/mov_INV1.gif'];
            system(processString);
            fprintf(fid,['INV1 image movie before/after correction:<br>\n']);
            fprintf(fid,['<img src="mov_INV1.gif"><br><br>\n']);
        end
        
        ov1 = orthoview(mOut.all_ims(:,:,:,2),'drawIms',0);
        imab_overwrite([htmlDir '/INV2.png'],ov1.oneIm);
        ov1 = orthoview(mOut.all_ims_corrected(:,:,:,2),'drawIms',0);
        imab_overwrite([htmlDir '/INV2_corrected.png'],ov1.oneIm);
        ov1 = orthoview(mOut.all_uniImage,'drawIms',0);
        imab_overwrite([htmlDir '/UNI.png'],ov1.oneIm);
        ov1 = orthoview(mOut.all_uniImage_corrected,'drawIms',0);
        imab_overwrite([htmlDir '/UNI_corrected.png'],ov1.oneIm);
        
        fprintf(fid,['INV2 image before correction:<br>\n']);
        fprintf(fid,['<img src="INV2.png"><br><br>\n']);
        fprintf(fid,['INV2 image after correction:<br>\n']);
        fprintf(fid,['<img src="INV2_corrected.png"><br><br>\n']);
        
        
        if testMagick==0
            processString = ['convert -dispose 2 -delay 50 -loop 0 ' htmlDir '/INV2.png ' htmlDir '/INV2_corrected.png ' htmlDir '/mov_INV2.gif'];
            system(processString);
            processString = ['convert -dispose 2 -delay 50 -loop 0 ' htmlDir '/UNI.png ' htmlDir '/UNI_corrected.png ' htmlDir '/mov_UNI.gif'];
            system(processString);
            fprintf(fid,['INV2 image movie before/after correction:<br>\n']);
            fprintf(fid,['<img src="mov_INV2.gif"><br><br>\n']);
            
        end
        
        fprintf(fid,['UNI image before correction:<br>\n']);
        fprintf(fid,['<img src="UNI.png"><br><br>\n']);
        fprintf(fid,['UNI image after correction:<br>\n']);
        fprintf(fid,['<img src="UNI_corrected.png"><br><br>\n']);
        
        if testMagick==0
            fprintf(fid,['UNI image movie before/after correction:<br>\n']);
            fprintf(fid,['<img src="mov_UNI.gif"><br><br>\n']);
        end
        
        % Do skull stripping with BET to be able to see the vessels within the brain more
        % clearly in the INV2 image
        disp('...............')
        disp('... Checking for FSL...')
        testFSL = getenv('FSLDIR');
        if ~isempty(testFSL)
            disp('... Found FSL, assuming that ''bet'' and ''fslmaths'' commands will work')
            disp('... Performing BET brain extraction ...')
            setenv('FSLOUTPUTTYPE','NIFTI'); % this uses more space than NIFTI_GZ, but stays compatible with SPM
            system(['bet ' outDir '/INV2_corrected ' outDir '/INV2_corrected_bet -m -f 0.2']);
            system(['fslmaths ' outDir '/INV2 -mul ' outDir '/INV2_corrected_bet_mask ' outDir '/INV2_bet']);
            disp('... Done')
            
            %%% make MIPs of INV2 after BET
            
            inv2 = rn([outDir '/INV2_bet.nii']);
            inv2c = rn([outDir '/INV2_corrected_bet.nii']);
            
            ov1 = orthoview(inv2,'mip',1,'drawIms',0);
            ov2 = orthoview(inv2c,'mip',1,'drawIms',0);
            
            im1 = abs(ov1.oneIm); im1 = im1/max(im1(:));
            im2 = abs(ov2.oneIm); im2 = im2/max(im2(:));
            
            imab_overwrite([htmlDir '/INV2_MIP.png'],im1);
            imab_overwrite([htmlDir '/INV2_corrected_MIP.png'],im2);
            
            if testMagick==0
                processString = ['convert -dispose 2 -delay 50 -loop 0 ' htmlDir '/INV2_MIP.png ' htmlDir '/INV2_corrected_MIP.png ' htmlDir '/mov_INV2_MIP.gif'];
                system(processString);
                
                fprintf(fid,['INV2 MIP movie before/after correction:<br>\n']);
                fprintf(fid,['<img src="mov_INV2_MIP.gif"><br><br>\n']);
            else
                fprintf(fid,['INV2 MIP before correction:<br>\n']);
                fprintf(fid,['<img src="INV2_MIP.png"><br><br>\n']);
                fprintf(fid,['INV2 MIP after correction:<br>\n']);
                fprintf(fid,['<img src="INV2_corrected_MIP.png"><br><br>\n']);
            end
            
        end
        
        %%% Add a zoom of the 3 images before and after correction
        if ~isempty(testFSL) % there will be a BET mask to use for choosing the brain region
            brainmask = rn([outDir '/INV2_corrected_bet_mask.nii']);            
            xi = round(Hxyz(1)/2+Hxyz(1)/10);
            yi = round(Hxyz(2)/2+Hxyz(2)/8:find(squeeze(any(any(brainmask,1),3)),1,'last'));
            zi = round(Hxyz(3)/2:find(squeeze(any(any(brainmask,1),2)),1,'last'));            
        else
            xi = round(Hxyz(1)/2+Hxyz(1)/10);
            yi = round(Hxyz(2)/2+Hxyz(2)/8: .9*Hxyz(2)); % arbitrary cut off points for zoom in y and z
            zi = round(Hxyz(3)/2:Hxyz(3): .8*Hxyz(3)); 
        end
        clim1 = percentile(mOut.all_ims(:,:,:,1),97);
        clim2 = percentile(mOut.all_ims(:,:,:,2),97);
        clim1c = percentile(mOut.all_ims_corrected(:,:,:,1),97);
        clim2c = percentile(mOut.all_ims_corrected(:,:,:,2),97);       
        clims_uni = [-.5 .5];
        
        fig(figIndex)
        clf
        set(gcf,'Position',[   246   611   982   494])
        subplot1(1,3)
        subplot1(1)
        imab(squeeze(mOut.all_ims(xi,yi,zi,1)),[0 clim1])
        subplot1(2)
        imab(squeeze(mOut.all_ims(xi,yi,zi,2)),[0 clim2])
        subplot1(3)
        imab(squeeze(mOut.all_uniImage(xi,yi,zi)),clims_uni)
        colormap(gray)
        export_fig([htmlDir '/zoom.png'])
        subplot1(1)
        imab(squeeze(mOut.all_ims_corrected(xi,yi,zi,1)),[0 clim1c])
        subplot1(2)
        imab(squeeze(mOut.all_ims_corrected(xi,yi,zi,2)),[0 clim2c])
        subplot1(3)
        imab(squeeze(mOut.all_uniImage_corrected(xi,yi,zi)),clims_uni)
        export_fig([htmlDir '/zoom_corrected.png'])
        if testMagick==0
            processString = ['convert -dispose 2 -delay 50 -loop 0 ' htmlDir '/zoom.png ' htmlDir '/zoom_corrected.png ' htmlDir '/mov_zoom.gif'];
            system(processString);
            
            fprintf(fid,['Sagittal zoom before/after correction:<br>\n']);
            fprintf(fid,['<img src="mov_zoom.gif"><br><br>\n']);
        else
            fprintf(fid,['Sagittal zoom before correction:<br>\n']);
            fprintf(fid,['<img src="zoom.png"><br><br>\n']);
            fprintf(fid,['Sagittal after correction:<br>\n']);
            fprintf(fid,['<img src="zoom_corrected.png"><br><br>\n']);
        end
        
        
        
        
end




%% And put a timing report into the html

stopTime = clock;
totalTime = etime(stopTime,startTime)/60/60;
totalTime_hrs = floor(totalTime);
if totalTime_hrs > 0
    totalTime_mins = round(rem(totalTime,totalTime_hrs)*60);
else
    totalTime_mins = round(totalTime*60);
end

fprintf(fid,['<h4>Total reconstruction time: ' num2str(totalTime_hrs) ' hours, ' num2str(totalTime_mins) ' mins</h4>\n']);
fprintf(fid,['<strong>Calculate GRAPPA weights for FatNavs: </strong>' num2str(round(timingReport_FatNavs.calculateGRAPPAweights)) ' seconds.<br>\n']);
fprintf(fid,['<strong>Reconstruct FatNavs: </strong>' num2str(nFatNavs) 'x ' num2str(round(mean(timingReport_FatNavs.eachFatNav))) ' seconds. Total time (possibly parallelized!): = ' num2str(round(timingReport_FatNavs.allFatNavs)) ' seconds. <br>\n']);
fprintf(fid,['<strong>Align FatNavs using SPM: </strong>' num2str(round(timingReport_FatNavs.SPMalignment)) ' seconds.<br>\n']);
if reconPars.bGRAPPAinRAM
    fprintf(fid,'<em>GRAPPA recon of host performed in RAM (i.e. faster)...</em><br>\n');
else
    fprintf(fid,'<em>GRAPPA recon of host stored as temporary files on hard drive (i.e. slower)...</em><br>\n');
end
if useGRAPPAforHost
    fprintf(fid,['<strong>Do 1D FFT for each ''slice'' in readout direction of host data:</strong> ' num2str(round(timingReport_hostRecon.FFTperSlice)) ' seconds.<br>\n']);
    fprintf(fid,['<strong>Declare variables for GRAPPA recon:</strong> ' num2str(round(timingReport_hostRecon.declareVariables)) ' seconds.<br>\n']);
    fprintf(fid,['<strong>GRAPPA recon:</strong> ' num2str(timingReport_hostRecon.GRAPPArecon/60,'%.1f') ' mins.<br>\n']);
end
fprintf(fid,['<strong>Application of retrospective motion-correction: </strong>' num2str(nc_keep) ' channels, ' num2str(nS) ' sets, each taking ' num2str(round(avgTimeApplyMocoPerVolume)) ' seconds (possibly parallelized) = ' num2str(timingReport_totalTimeApplyMoco/60,'%.1f') ' mins.<br>\n']);

% include version number
fprintf(fid,['<br><br><br><em>' char(datetime) '- created with reconstructSiemensMP2RAGEwithFatNavs.m, version: ' retroMocoBoxVersion '</em>\n']);


fprintf(fid,'</body></html>\n');
fclose(fid);


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

if ~reconPars.bKeepFatNavs
    rmdir(fatnavdir,'s')
end
    
timingReport.totalTime_hrs = totalTime_hrs;
timingReport.totalTime_mins = totalTime_mins;
timingReport.timingReport_FatNavs = timingReport_FatNavs;
timingReport.timingReport_hostRecon = timingReport_hostRecon;
timingReport.avgTimeApplyMocoPerVolume = avgTimeApplyMocoPerVolume;
timingReport.timingReport_totalTimeApplyMoco = timingReport_totalTimeApplyMoco;


end

