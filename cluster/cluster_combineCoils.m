% function cluster_combineCoils(reconParsFile)

load(reconParsFile)

run([RETROMOCOBOX_PATH '/addRetroMoCoBoxToPath.m']);
% run([MIRT_PATH '/setup.m']);
addpath(SPM_PATH);

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


all_ims = zeros(Hxyz(1),Hxyz(2),Hxyz(3),nEco,'single'); % assume this *will* fit in RAM - ends up at 1.2 GB for 336x336x192x7
all_ims_corrected = zeros(Hxyz(1),Hxyz(2),Hxyz(3),nEco,'single'); % two copies required if we want to keep uncorrected as well...!

if nEco > 1
    phase_diff = complex(zeros(Hxyz(1),Hxyz(2),Hxyz(3),'single')); 
    phase_diff_corrected = complex(zeros(Hxyz(1),Hxyz(2),Hxyz(3),'single')); 
end

parfor iC = 1:nc_keep
    
    thisCoil = load([tempNameRoots.finalImage '_' num2str(iC,'%.2d') '.mat']);
    
    all_ims = all_ims + abs(single(thisCoil.thiscoil_mag)).^2;
    all_ims_corrected = all_ims_corrected + abs(single(thisCoil.thiscoil_mag_corrected)).^2;
    
    if nEco > 1
        phase_diff = phase_diff + abs(single(thisCoil.thiscoil_mag(:,:,:,1))).*exp(1i*2*pi*((single(thisCoil.thisphase_diff)/4095)-.5));
        phase_diff_corrected = phase_diff_corrected + abs(single(thisCoil.thiscoil_mag_corrected(:,:,:,1))).*exp(1i*2*pi*((single(thisCoil.thisphase_diff_corrected)/4095)-.5));
    end 
end

                      
all_ims = sqrt(all_ims);
all_ims_corrected = sqrt(all_ims_corrected);

sn(int16(4095*all_ims/max(reshape(all_ims,[],1))),[outDir '/GRE'],hostVoxDim_mm)
sn(int16(4095*all_ims_corrected/max(reshape(all_ims_corrected,[],1))),[outDir '/GRE_corrected'],hostVoxDim_mm)
    
if nEco > 1
    sn(int16(4095* (angle(phase_diff)+pi)/(2*pi)),[outDir '/GRE_phasediff'],hostVoxDim_mm)
    sn(int16(4095* (angle(phase_diff_corrected)+pi)/(2*pi)),[outDir '/GRE_phasediff_corrected'],hostVoxDim_mm)
end

