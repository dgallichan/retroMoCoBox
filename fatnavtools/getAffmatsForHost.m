function [this_fitMat_mm, thisRot, rotMatDisplacement_mm] = getAffmatsForHost(this_fitMat,twix_obj,reconPars)

% this_fitMat = fitResult.MPos_cent.mats;
% reconPars.swapDims_xyz

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

%%






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




% if reconPars.swapDims_xyz(alignDim)
%     % make it so that centre of k-space is not 'moved' (accounting for partial Fourier):
%     this_fitMat_mm = recentre_affmats(this_fitMat_mm,hxyz(alignDim)-kspaceCentre_xyz(alignDim));  
%     % alignIndices = hxyz(alignDim):-1:1;  
% else
%     % make it so that centre of k-space is not 'moved' (accounting for partial Fourier):
%     this_fitMat_mm = recentre_affmats(this_fitMat_mm,kspaceCentre_xyz(alignDim));  
%     % alignIndices = 1:hxyz(alignDim);
% end