function cluster_runMultiEchoGRE_eachCoil(reconParsFile,iC,bSaveAll)

if nargin < 3
    bSaveAll = 0;
end

load(reconParsFile)

run([RETROMOCOBOX_PATH '/addRetroMoCoBoxToPath.m']);
% run([MIRT_PATH '/setup.m']);
addpath(SPM_PATH);

% if nEco <= 10
%     parpool(nEco); % this line fails on the cluster for some reason...!
% end

thiscoil_ims = complex(zeros(Hxyz(1),Hxyz(2),Hxyz(3),nEco,'single')); % assume this *will* fit in RAM - ends up at 1.2 GB for 336x336x192x7
thiscoil_ims_corrected = complex(zeros(Hxyz(1),Hxyz(2),Hxyz(3),nEco,'single')); % two copies required if we want to keep uncorrected as well...!

% % reshape the array to see if parfor can slice it properly then...
% (doesn't seem to work at reducing RAM requirements for this script -
% still needs ~30 GB for 336x336x192x7 data, which makes no sense!)
% thiscoil_ims = complex(zeros(nEco,Hxyz(1),Hxyz(2),Hxyz(3),'single')); % assume this *will* fit in RAM - ends up at 1.2 GB for 336x336x192x7
% thiscoil_ims_corrected = complex(zeros(nEco,Hxyz(1),Hxyz(2),Hxyz(3),'single')); % two copies required if we want to keep uncorrected as well...!


% parfor doesn't like some stuff... have to make it clear that they are
% variables by setting them equal to themselves...!
tempGRAPPAreconFile = tempNameRoots.grappaRecon_1DFFT;
swapDims_xyz = reconPars.swapDims_xyz;
nc_keep = nc_keep;
hrps = hrps;
iS = iS;
iC_keep = iC_keep;
permutedims = permutedims;
hxyz = hxyz;
Hxyz = Hxyz;
nEco = nEco;
fitMats_mm_toApply = fitMats_mm_toApply;
alignDim = alignDim;
alignIndices = alignIndices;
hostVoxDim_mm = hostVoxDim_mm;
kspaceCentre_xyz = kspaceCentre_xyz;
tempNameRoots = tempNameRoots;
MIDstr = MIDstr;

parfor iE = 1:nEco
    
    fprintf(['Reconstructing coil ' num2str(iC) ' of ' num2str(nc_keep) ', echo ' num2str(iE) '\n']);
    
    
    thisData = zeros(hrps.');
    for iReadSlice = 1:hrps(1) % virtual 'slices' in the readout direction
        tempData = load([tempGRAPPAreconFile '_' num2str(iReadSlice) '_' num2str(iE) '_' num2str(iS) '.mat']);
        thisData(iReadSlice,:,:) = reshape(tempData.outData(1,iC_keep(iC),:,:),[1 hrps(2) hrps(3)]);
    end
    
    thisData = fft1s(thisData,1); % put into full 3D k-space
    
    
    thisData = permute(thisData,permutedims);
    if swapDims_xyz(1)
        thisData = thisData(end:-1:1,:,:,:);
    end
    if swapDims_xyz(2)
        thisData = thisData(:,end:-1:1,:,:);
    end
    if swapDims_xyz(3)
        thisData = thisData(:,:,end:-1:1,:);
    end
    
    newData = thisData;
    
    if any(hxyz~=Hxyz)
        thisData(Hxyz(1),Hxyz(2),Hxyz(3)) = 0; % extend to new size
        thisData = circshift(thisData,double([Hxyz(1)-hxyz(1) Hxyz(2)-hxyz(2) Hxyz(3)-hxyz(3)].*(1-swapDims_xyz)));
        % the 'double' in the above line appears necessary in certain
        % versions of Matlab, no idea why...
    end
    
    thisData = ifft3s(thisData)*prod(hxyz);
    
    thisData_corrected = applyRetroMC_nufft(newData,fitMats_mm_toApply,alignDim,alignIndices,11,hostVoxDim_mm,Hxyz,kspaceCentre_xyz);
    
    
    thiscoil_ims(:,:,:,iE) = thisData;
    thiscoil_ims_corrected(:,:,:,iE) = thisData_corrected;

%     thiscoil_ims(iE,:,:,:) = reshape(thisData,[1 Hxyz.']);
%     thiscoil_ims_corrected(iE,:,:,:) = reshape(thisData_corrected,[1 Hxyz.']);
    
end

% thiscoil_ims = permute(thiscoil_ims,[2 3 4 1]);
% thiscoil_ims_corrected = permute(thiscoil_ims_corrected,[2 3 4 1]);

if bSaveAll % this uses more disk space, but can be plugged in more easily to stuff like ASPIRE coil combination / phase reconstruction
    save([tempNameRoots.finalImage '_' num2str(iC,'%.2d') '.mat'],'thiscoil_ims','thiscoil_ims_corrected');
    
else

    % After motion-correction has been applied and we are working only in
    % image space we can afford to reduce the precision of the data, but it
    % needs to be normalised the same for each coil... It may not always be
    % true, but for the test data it looks like it is scaled on the order
    % of magnitude of one - so this should be OK
    thiscoil_mag = int16(4095*abs(thiscoil_ims));
    thiscoil_mag_corrected = int16(4095*abs(thiscoil_ims_corrected));
    
    if nEco > 1
        thisphase_diff = int16((4095/(2*pi))*(angle(mean(thiscoil_ims(:,:,:,2:end).*conj(thiscoil_ims(:,:,:,1:end-1)),4))+pi) );
        thisphase_diff_corrected =  int16((4095/(2*pi))*(angle(mean(thiscoil_ims_corrected(:,:,:,2:end).*conj(thiscoil_ims_corrected(:,:,:,1:end-1)),4))+pi));
        save([tempNameRoots.finalImage '_' num2str(iC,'%.2d') '.mat'],'thiscoil_mag','thiscoil_mag_corrected','thisphase_diff','thisphase_diff_corrected');
    else
        save([tempNameRoots.finalImage '_' num2str(iC,'%.2d') '.mat'],'thiscoil_mag','thiscoil_mag_corrected');
    end
end