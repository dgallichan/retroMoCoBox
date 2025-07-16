function cluster_runMultiEchoGRE_eachEcho_skipMoCo(reconParsFile,iE)

% _skipMoCo version strips out all MoCo to speedup getting through ASPIRE
% pipeline for testing

load(reconParsFile)

run([RETROMOCOBOX_PATH '/addRetroMoCoBoxToPath.m']);
% run([MIRT_PATH '/setup.m']);
addpath(SPM_PATH);

% if nEco <= 10
%     parpool(nEco); % this line fails on the cluster for some reason...!
% end
isPool = gcp('nocreate');
if isempty(isPool)
    c = parcluster('local');
    c.NumWorkers = reconPars.parpoolSize;
    mkdir(tempDir,['matlab_tmp_' num2str(iE)]);
    % use a specific job storate location to try to avoid problems with
    % multiple parpools simultaneously running on cluster
    c.JobStorageLocation = [tempDir '/matlab_tmp_' num2str(iE)];
    c.parpool(c.NumWorkers);
end


% % reshape the array to see if parfor can slice it properly then...
all_ims = complex(zeros(nc_keep,Hxyz(1),Hxyz(2),Hxyz(3),'single')); 
% all_ims_corrected = complex(zeros(nc_keep,Hxyz(1),Hxyz(2),Hxyz(3),'single')); % two copies required if we want to keep uncorrected as well...!


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

tempDir = fileparts(reconParsFile); % tempDir itself wasn't actually saved explicitly before...


% Setup the objects for the NUFFT:
% oversampFactor = 1.5;
% [~, st, ~, phaseTranslations] = applyRetroMC_nufft(zeros(hxyz'),fitMats_mm_toApply,alignDim,alignIndices,11,hostVoxDim_mm,Hxyz,kspaceCentre_xyz,-1, oversampFactor, 1);

% load ALL the data into RAM for each echo
allGRAPPAdata = zeros(nc_keep,hrps(1),hrps(2),hrps(3));
for iReadSlice = 1:hrps(1)% virtual 'slices' in the readout direction
    tempData = load([tempGRAPPAreconFile '_' num2str(iReadSlice) '_' num2str(iE) '_' num2str(iS) '.mat']);
    allGRAPPAdata(:,iReadSlice,:,:) = reshape(tempData.outData,[nc_keep 1 hrps(2) hrps(3)]);
end
        

parfor iC = 1:nc_keep
        
        fprintf(['Reconstructing coil ' num2str(iC) ' of ' num2str(nc_keep) ', echo ' num2str(iE) '\n']);
        
        
        thisData = squeeze(allGRAPPAdata(iC,:,:,:));
        
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
        
%         thisData_corrected = applyRetroMC_nufft(newData,fitMats_mm_toApply,alignDim,alignIndices,11,hostVoxDim_mm,Hxyz,kspaceCentre_xyz);
        % thisData_corrected = nufft_adj_single(newData(:).*phaseTranslations(:),st);        
        
        all_ims(iC,:,:,:) = reshape(thisData,[1 Hxyz.']);
        % all_ims_corrected(iC,:,:,:) = reshape(thisData_corrected,[1 Hxyz.']);


end

all_ims = permute(all_ims,[2 3 4 1]);
% all_ims_corrected = permute(all_ims_corrected,[2 3 4 1]);

% create 4D NIFTI files [x y z coils]
sn(abs(all_ims),[tempDir '/GRE_Echo' num2str(iE) '_mag'],hostVoxDim_mm) 
sn(angle(all_ims),[tempDir '/GRE_Echo' num2str(iE) '_phs'],hostVoxDim_mm)

% sn(abs(all_ims_corrected),[tempDir '/GRE_MoCo_Echo' num2str(iE) '_mag'],hostVoxDim_mm)
% sn(angle(all_ims_corrected),[tempDir '/GRE_MoCo_Echo' num2str(iE) '_phs'],hostVoxDim_mm)
% 
% % create 5D NIFTI files [x y z echoes coils]
% sn(abs(all_ims),[tempDir '/GRE_5D_mag'],hostVoxDim_mm) 
% sn(angle(all_ims),[tempDir '/GRE_5D_phs'],hostVoxDim_mm)
% 
% sn(abs(all_ims_corrected),[tempDir '/GRE_5D_MoCo_mag'],hostVoxDim_mm)
% sn(angle(all_ims_corrected),[tempDir '/GRE_5D_MoCo_phs'],hostVoxDim_mm)
