% function cluster_combineCoils(reconParsFile)

load(reconParsFile)

run([RETROMOCOBOX_PATH '/addRetroMoCoBoxToPath.m']);
% run([MIRT_PATH '/setup.m']);
addpath(SPM_PATH);

tempDir = fileparts(reconParsFile); % tempDir itself wasn't actually saved explicitly before...

all_ims = complex(zeros(Hxyz(1),Hxyz(2),Hxyz(3),nEco,nc_keep,'single')); % assume this *will* fit in RAM - ends up at 38 GB for 336x336x192x7x32
all_ims_corrected = complex(zeros(Hxyz(1),Hxyz(2),Hxyz(3),nEco,nc_keep,'single')); % two copies required if we want to keep uncorrected as well...!


disp('Loading all echoes into RAM...')
for iE = 1:nEco
    
    thisEco = rn([tempDir '/GRE_Echo' num2str(iE) '_mag']).*exp(1i*rn([tempDir '/GRE_Echo' num2str(iE) '_phs']));
    thisEco_corrected = rn([tempDir '/GRE_MoCo_Echo' num2str(iE) '_mag']).*exp(1i*rn([tempDir '/GRE_MoCo_Echo' num2str(iE) '_phs']));
    
    all_ims(:,:,:,iE,:) = reshape(thisEco,[Hxyz' 1 nc_keep]);
    all_ims_corrected(:,:,:,iE,:) = reshape(thisEco_corrected,[Hxyz' 1 nc_keep]); 
    
    disp(iE)
end
disp('Done!')


disp('Saving out 4 combined NIFTI files...')
% create 5D NIFTI files [x y z echoes coils]
sn(abs(all_ims),[tempDir '/GRE_5D_mag'],hostVoxDim_mm) 
disp('1')
sn(angle(all_ims),[tempDir '/GRE_5D_phs'],hostVoxDim_mm)
disp('2')

sn(abs(all_ims_corrected),[tempDir '/GRE_5D_MoCo_mag'],hostVoxDim_mm)
disp('3')
sn(angle(all_ims_corrected),[tempDir '/GRE_5D_MoCo_phs'],hostVoxDim_mm)
disp('4')
disp('Done!')

