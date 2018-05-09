function reconstructSiemensGREwithFatNavs_cluster(rawDataFile, varargin)

[reconPars.outRoot, reconPars.tempRoot, reconPars.bLinParSwap, reconPars.bGRAPPAinRAM, reconPars.bKeepGRAPPArecon, reconPars.bKeepReconInRAM, reconPars.bFullParforRecon,...
    reconPars.coilCombineMethod, reconPars.FatNavRes_mm, reconPars.swapDims_xyz, reconPars.bZipNIFTIs, reconPars.bKeepFatNavs,reconPars.bKeepPatientInfo,...
    CLUSTER_LOG_PATH] = process_options(varargin,...
    'outRoot',[],'tempRoot',[],'bLinParSwap',0,'bGRAPPAinRAM',0,'bKeepGRAPPArecon',0,'bKeepReconInRAM',0,...
    'bFullParforRecon',0,'coilCombineMethod','default','FatNavRes_mm',[],'swapDims_xyz',[0 0 1],'bZipNIFTIs',1,'bKeepFatNavs',0,'bKeepPatientInfo',1,...
    'CLUSTER_LOG_PATH','~/');

% bKeepReconInRAM, bFullParforRecon, coilCombineMethod - all ignored for this GRE script!!

reconPars.retroMocoBoxVersion = '0.7.0dev';



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
twix_obj = mapVBVD_fatnavs(rawDataFile,'removeOS',1);
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
reconPars.FatNavRes_mm = twix_obj.hdr.MeasYaps.sWiPMemBlock.alFree{5};



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

[datadir, fileName, fileExt] = fileparts(rawDataFile);
if isempty(reconPars.outRoot)
    reconPars.outRoot = datadir;
end
if reconPars.outRoot(1)=='~'  % SPM seems to have a 'thing' about the tilde...
    reconPars.outRoot = [getenv('HOME') reconPars.outRoot(2:end)];
end

figIndex = 999; % figure to use for outputs

MIDstr = getMIDstr(rawDataFile);

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
fidHTML = fopen([htmlDir '/index.html'],'w');
fprintf(fidHTML,['<html><head><title>GRE with FatNavs - Summary</title>\n']);
fprintf(fidHTML,'</head>\n');
fprintf(fidHTML,'<body>\n');
fprintf(fidHTML,['<h2>GRE with FatNavs - Summary: %s</h2>\n'],[fileName]);
fprintf(fidHTML,['<br>\n']);





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
nEco = twix_obj.image.NEco;

fprintf(fidHTML,['<h4>Host GRE sequence</h4>\n']);
fprintf(['\n\n\nHost GRE sequence\n']);

fprintf(fidHTML,['RPS dimensions in raw data: ' num2str(hrps(1)) 'x' num2str(hrps(2)) 'x' num2str(hrps(3)) '<br>\n']);
fprintf(['RPS dimensions in raw data: ' num2str(hrps(1)) 'x' num2str(hrps(2)) 'x' num2str(hrps(3)) '\n']);
fprintf(fidHTML,['XYZ dimensions reconstructed: ' num2str(Hxyz(1)) 'x' num2str(Hxyz(2)) 'x' num2str(Hxyz(3)) '<br>\n']);
fprintf(['XYZ dimensions reconstructed: ' num2str(Hxyz(1)) 'x' num2str(Hxyz(2)) 'x' num2str(Hxyz(3)) '\n']);
fprintf(fidHTML,['FOV - ' num2str(FOVxyz(1),'%.1f') 'x' num2str(FOVxyz(2),'%.1f') 'x' num2str(FOVxyz(3),'%.1f') 'mm<br>\n']);
fprintf(['FOV - ' num2str(FOVxyz(1),'%.1f') 'x' num2str(FOVxyz(2),'%.1f') 'x' num2str(FOVxyz(3),'%.1f') 'mm\n']);
fprintf(fidHTML,['Resolution: ' num2str(hostVoxDim_mm(1),'%.3f') 'x' num2str(hostVoxDim_mm(2),'%.3f') 'x' num2str(hostVoxDim_mm(3),'%.3f') 'mm<br>\n']);
fprintf(['Resolution: ' num2str(hostVoxDim_mm(1),'%.3f') 'x' num2str(hostVoxDim_mm(2),'%.3f') 'x' num2str(hostVoxDim_mm(3),'%.3f') 'mm\n']);
if manualFatNavRes
    fprintf(fidHTML,['Manually selected FatNav resolution: ' num2str(reconPars.FatNavRes_mm) ' mm<br>\n']);
    fprintf(['Manually selected FatNav resolution: ' num2str(reconPars.FatNavRes_mm) ' mm\n']);
else
    fprintf(fidHTML,['Assumed FatNav resolution (based on field strength): ' num2str(reconPars.FatNavRes_mm) ' mm<br>\n']);
    fprintf(['Assumed FatNav resolution (based on field strength): ' num2str(reconPars.FatNavRes_mm) ' mm\n']);
end
fprintf(fidHTML,['Detected orientation: ' orientText '<br>\n']);
fprintf(['Detected orientation: ' orientText '\n']);
if isfield(twix_obj.hdr.MeasYaps,'sCoilSelectMeas')
    fprintf(['Coil used: ' char(twix_obj.hdr.MeasYaps.sCoilSelectMeas.aRxCoilSelectData{1}.asList{1}.sCoilElementID.tCoilID) ', with ' num2str(nc) ' channels active\n']);
    fprintf(fidHTML,['Coil used: ' char(twix_obj.hdr.MeasYaps.sCoilSelectMeas.aRxCoilSelectData{1}.asList{1}.sCoilElementID.tCoilID) ', with ' num2str(nc) ' channels active<br>\n']);
else
    fprintf(['Coil used: ' char(twix_obj.hdr.MeasYaps.asCoilSelectMeas{1}.asList{1}.sCoilElementID.tCoilID) ', with ' num2str(nc) ' channels active\n']);
    fprintf(fidHTML,['Coil used: ' char(twix_obj.hdr.MeasYaps.asCoilSelectMeas{1}.asList{1}.sCoilElementID.tCoilID) ', with ' num2str(nc) ' channels active<br>\n']);
end

%% Check if using HEADNECK_64 receive coil, and discard channels over the neck if this is the case (would be nice to know what Siemens does...)
if isfield(twix_obj.hdr.MeasYaps,'sCoilSelectMeas') ...
        && strcmp(twix_obj.hdr.MeasYaps.sCoilSelectMeas.aRxCoilSelectData{1}.asList{1}.sCoilElementID.tCoilID,'"HeadNeck_64"')
    iC_keep = 1:nc;
    iC_keep([1 7 8 18 29 30 39 40 49 50]) = []; % these channels cover the neck (in one test dataset with 52 data channels from the 64-channel coil...)
    fprintf(['\n\n****  Detected use of HeadNeck_64 RF coil **** \n'...
                 'Using manually predefined set of channels\n'...
                 'to reduce signal from neck area\n' ...
                 '*********************************\n\n']);
    fprintf(fidHTML,['<br><br>\n\n****  Detected use of HeadNeck_64 RF coil **** <br>\n'...
                 'Using manually predefined set of channels<br>\n'...
                 'to reduce signal from neck area<br>\n' ...
                 '*********************************<br><br>\n\n']);    
else
    iC_keep = 1:nc;
end

nc_keep = length(iC_keep);
    



%% Check the number of FatNavs available compared to the size of the host data

nFatNavs = twix_obj.FatNav.dataSize(9)/twix_obj.image.NAve;

alignDim = permutedims_toXYZ(2);
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
            reconPars.outRoot,'FatNavRes_mm',reconPars.FatNavRes_mm, 'iAve', reconPars.iAve, 'appendString', appendString);

% And put stuff into the html report
if exist(fatnavdir,'dir') % could have just kept the motion-parameters file...
    imdims = 2*[304 128];
    
    if exist([fatnavdir '/a_FatNav_ACSim_' MIDstr '.png'],'file')
        copyfile([fatnavdir '/a_FatNav_ACSim_' MIDstr '.png'],[htmlDir '/ACSim.png']);
        fprintf(fidHTML,['ACS image:<br>\n']);
        fprintf(fidHTML,['<img src="ACSim.png" height=%s width=%s><br><br>\n'],num2str(imdims(2)),num2str(imdims(1)));
    end
    if exist([fatnavdir '/a_FatNav_MoCoPars_' MIDstr '.png'],'file')
        copyfile([fatnavdir '/a_FatNav_MoCoPars_' MIDstr '.png'],[htmlDir '/motion_parameters.png']);
        fprintf(fidHTML,['<h4>FatNavs</h4>\n']);
        fprintf(fidHTML,['Estimated motion parameters:<br>\n']);
        fprintf(fidHTML,['<img src="motion_parameters.png"><br><br>\n']);
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
title(['xy MIP of host GRAPPA recon, coil ' num2str(iAsymCoil)])
colormap(gray)
export_fig([htmlDir '/orientationCheck_xy.png']);
fprintf(fidHTML,['Orientation check for left/right symmetry:<br>\n']);
fprintf(fidHTML,['<img src="orientationCheck_xy.png"><br><br>\n']);
fprintf(fidHTML,['(Both images above should have the brightest signal on the left of the image. If not, the orientation of the FatNavs is not correctly aligned with the host sequence)<br><br><br>\n']);


%%

if exist([fatnavdir '/eachFatNav_001.nii'],'file')
    affMat = thisRot;
    affMat(1:3,4) = rotMatDisplacement_mm.'/2;
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
    fprintf(fidHTML,['Orientation check for host sequence slice rotation and positioning:<br>\n']);
    fprintf(fidHTML,['<img src="orientationCheck_FatVolume.png"><br><br>\n']);
    fprintf(fidHTML,['(The fat volume shown above should approximately correspond to the FOV chosen for the host sequence)<br><br><br>\n']);
end

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



%% Apply the motion-correction by spawning jobs on the cluster
% (also need to spawn clean-up and finishing jobs which wait for this to
% finish)

RETROMOCOBOX_PATH = fileparts(which('addRetroMoCoBoxToPath.m'));
mirtpath = fileparts(which('nufft.m'));
MIRT_PATH = [mirtpath '/../'];
SPM_PATH = fileparts(which('spm.m')); 

tempNameRoots.allReconPars = [tempDir '/tempAllReconPars_' MIDstr '.mat'];
tempNameRoots.finalImage = [tempDir '/tempFinalImage_' MIDstr];
tempNameRoots.clusterScript = [tempDir '/tempClusterScript_' MIDstr '.sh'];
tempNameRoots.clusterScriptRecombine = [tempDir '/tempClusterScriptRecombine_' MIDstr '.sh'];
tempNameRoots.clusterScriptCleanup = [tempDir '/tempClusterScriptCleanup_' MIDstr '.sh'];
tempNameRoots.cleanupFiles = [tempDir '/tempCleanupPars_' MIDstr '.mat'];

iS = 1; % GRE data has only one 'set'
save(tempNameRoots.allReconPars,'iS','nc_keep','Hxyz','hxyz','nEco','reconPars','iC_keep','hrps','permutedims',...
'fitMats_mm_toApply','alignDim','alignIndices','hostVoxDim_mm','kspaceCentre_xyz','tempNameRoots','MIDstr','outDir',...
'RETROMOCOBOX_PATH','MIRT_PATH','SPM_PATH')

fid = fopen(tempNameRoots.clusterScript,'w');
fprintf(fid,'#!/bin/bash\n');
fprintf(fid,['#SBATCH --array 1-' num2str(nc) '\n']);
fprintf(fid,'#SBATCH -p cubric-default\n');
fprintf(fid,'#SBATCH --job-name=GREreconHelper\n');
fprintf(fid,['#SBATCH -o ' CLUSTER_LOG_PATH '/GREreconHelperArray_%j.out\n']);
fprintf(fid,['#SBATCH -e ' CLUSTER_LOG_PATH '/GREreconHelperArray_%j.err\n']);
if nEco <= 10
    fprintf(fid,['#SBATCH --ntasks ' num2str(nEco) '\n']); % only 7 echoes, so parfor only needs 7 threads
else
    fprintf(fid,'#SBATCH --ntasks 10\n'); 
end
fprintf(fid,'#SBATCH --mem-per-cpu=32000M\n'); % MaxRSS showed requiring ~30 GB RAM for 336x336x192x7 data - I don't know why it still needs so much...
%%%% I even tried switching around the parfor loop in cluster_runMultieEchoGRE_eachcoil to try to get the RAM usage down - I'm not sure if this is a bug in this version of MATLAB...
fprintf(fid,['cd ' RETROMOCOBOX_PATH '/cluster\n']);
fprintf(fid,['matlab -nodesktop -nosplash -r "cluster_runMultiEchoGRE_eachCoil(''' tempNameRoots.allReconPars ''',${SLURM_ARRAY_TASK_ID});exit;"\n']);
fclose(fid);

disp('Launching batch job array for applying the motion correction...')
[status, sbatch_out] = system(['sbatch ' tempNameRoots.clusterScript]); 
jobnumber = sbatch_out(21:end); % after the text 'Submitted batch job ..'
disp('Done.')

fid = fopen(tempNameRoots.clusterScriptRecombine,'w');
fprintf(fid,'#!/bin/bash\n');
fprintf(fid,['#SBATCH --dependency=afterok:' jobnumber]); % wait for the array above to finish before starting the recombine
fprintf(fid,'#SBATCH -p cubric-default\n');
fprintf(fid,'#SBATCH --job-name=GREreconRecombine\n');
fprintf(fid,['#SBATCH -o ' CLUSTER_LOG_PATH '/GREreconRecombine_%j.out\n']);
fprintf(fid,['#SBATCH -e ' CLUSTER_LOG_PATH '/GREreconRecombine_%j.err\n']);
fprintf(fid,'#SBATCH --ntasks 10\n');
fprintf(fid,'#SBATCH --mem-per-cpu=32000M\n'); % this one also exceeded memory limit when set to 16000
fprintf(fid,['cd ' RETROMOCOBOX_PATH '/cluster\n']);
fprintf(fid,['matlab -nodesktop -nosplash -r "reconParsFile = ''' tempNameRoots.allReconPars ''';cluster_combineCoils;exit;"\n']);
fclose(fid);

disp('Launching batch job array for applying the motion correction...')
[status, sbatch_out] = system(['sbatch ' tempNameRoots.clusterScriptRecombine]); 
jobnumber2 = sbatch_out(21:end); % after the text 'Submitted batch job ..'
disp('Done.')

save(tempNameRoots.cleanupFiles,'htmlDir','startTime','timingReport_FatNavs','nc_keep','nEco', ...
    'reconPars','tempNameRoots','tempDir','outDir','nFatNavs','fatnavdir','RETROMOCOBOX_PATH','MIRT_PATH','SPM_PATH');
fclose(fidHTML); % HTML needs closing and will be reopened in cleanup

fid = fopen(tempNameRoots.clusterScriptCleanup,'w');
fprintf(fid,'#!/bin/bash\n');
fprintf(fid,['#SBATCH --dependency=afterok:' jobnumber2]); % wait for the recombine to finish before starting the cleanup
fprintf(fid,'#SBATCH -p cubric-default\n');
fprintf(fid,'#SBATCH --job-name=GREreconCleanup\n');
fprintf(fid,['#SBATCH -o ' CLUSTER_LOG_PATH '/GREreconCleanup_%j.out\n']);
fprintf(fid,['#SBATCH -e ' CLUSTER_LOG_PATH '/GREreconCleanup_%j.err\n']);
fprintf(fid,'#SBATCH --ntasks 1\n');
fprintf(fid,'#SBATCH --mem-per-cpu=8000M\n');
fprintf(fid,['cd ' RETROMOCOBOX_PATH '/cluster\n']);
fprintf(fid,['matlab -nodesktop -nosplash -r "cleanupFile = ''' tempNameRoots.cleanupFiles ''';cluster_cleanup;exit;"\n']);
fclose(fid);

disp('Launching batch job array for applying the motion correction...')
[status, sbatch_out] = system(['sbatch ' tempNameRoots.clusterScriptCleanup]); 
disp('Done.')







