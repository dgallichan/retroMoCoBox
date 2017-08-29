% run_retroMoCoDemo_simulateMotion.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Whereas run_retroMoCoDemo.m uses real example data acquired with a pulse
% sequence with FatNavs inserted into it to demonstrate the
% motion-correction, this script uses the retroMoCoBox to simulate
% motion-corrupted data - and then see how well the corrupted k-space data
% can be fixed by applying the known translations and rotations.
%
% This can't be considered a 'full' simulation of the artifacts caused by
% motion - as it makes a number of simplifying assumptions - but it should
% give a good indication of the upper limit of the quality of a
% retrospective motion-correction when the motion-parameters are known
% precisely (which for any real acquisition they will not be - and for a
% FatNav acquisition the motion parameters will necessarily not correspond
% exactly in time to when the k-space data for the main sequence were
% acquired).
%
% It's probably easiest to run this script cell-by-cell so that you can
% keep track of what's going on. The default settings should demonstrate
% that even for quite considerable motion, the single-step retrospective
% correction can do a good job if the motion-parameters are known exactly.
%
% The script lets you experiment with various types of motion profiles -
% with simple models for swallowing artifacts and sudden jerky movements.
% You can also vary the 'baseline' noise. If you change this from the
% 'smoother' motion (default) to the '*really* rough motion' you should
% find that although the motion-corrected version is reasonably sharp,
% there is still quite a lot of residual ringing artifact. This residual
% ringing can be reduced considerably by an iterative version of the
% reconstruction - but this is *much* slower, so should only be attempted
% on a reasonable small test volume (easiest to achieve by uncommenting the
% line in the next cell which will extract just a thin slab in the z
% direction). The iterative reconstruction massively reduces the ringing
% artifact, but still leaves a small textural artifact. It is an ongoing
% area of research as to what extent the iterative reconstruction can be
% used on real data.
%
%
% -- Daniel Gallichan, gallichand@cardiff.ac.uk, August 2017

%% -- SET PATHS MANUALLY IN THIS SECTION -- 
run('addRetroMoCoBoxToPath.m')

% The NUFFT uses the Michigan Image Reconstruction Toolbox (MIRT)
% (http://web.eecs.umich.edu/~fessler/code/index.html)
run('../mirt/setup.m')

%%% Load in an example image: 
%%% (The Colin27 brain is good for this - downloadable from here: http://www.bic.mni.mcgill.ca/ServicesAtlases/Colin27) 
image_original = rn('../exampleData/colin27_t1_tal_lin.nii'); hostVoxDim_mm = [1 1 1]; 
% image_original = rn('/usr/local/fsl/data/standard/MNI152_T1_2mm.nii.gz'); hostVoxDim_mm = [2 2 2];

% force dimension to be even for simplicity of consitent indexing:
[nx,ny,nz] = size(image_original);
nx = 2*floor(nx/2); ny = 2*floor(ny/2); nz = 2*floor(nz/2);
image_original = image_original(1:nx,1:ny,1:nz);

% image_original = image_original(:,:,80:100); % <--- use only a subset of the data to be much faster

% normalize:
image_original = image_original / percentile(abs(image_original),95);

rawData = fft3s(image_original);

nT = size(rawData,2);

%% Generate the artificial motion parameters - the magnitude of different components can be varied by hand!

rng(1); % Set the seed for the random number generator to be able to create reproducible motion patterns

% for Perlin noise, this determines the weights between different harmonics of noise
% noiseBasePars = 1; %% *really* rough motion
% noiseBasePars = 5;  %% quite 'rough' motion
noiseBasePars = 3.^[0:8]; %% smoother motion

maxDisp = 4; % magnitude of general background noise movement - translations
maxRot = 4; % magnitude of rotations
swallowFrequency = 3; % number of swallowing events in scan
swallowMagnitude = [3 3]; % first is translations, second is rotations 
suddenFrequency = 5; % number of sudden movements
suddenMagnitude = [3 3]; % first is translations, second is rotations 


% general background noise movement:
fitpars = maxDisp*(perlinNoise1D(nT,noiseBasePars).'-.5);
fitpars(2,:) = maxDisp*(perlinNoise1D(nT,noiseBasePars).'-.5);
fitpars(3,:) = maxDisp*(perlinNoise1D(nT,noiseBasePars).'-.5);

fitpars(4,:) = maxRot*(perlinNoise1D(nT,noiseBasePars).'-.5);
fitpars(5,:) = maxRot*(perlinNoise1D(nT,noiseBasePars).'-.5);
fitpars(6,:) = maxRot*(perlinNoise1D(nT,noiseBasePars).'-.5);

% add in swallowing-like movements - just to z direction and pitch:
swallowTraceBase = exp(-linspace(0,1e2,nT));
swallowTrace = zeros(1,nT);
for iS = 1:swallowFrequency
    swallowTrace = swallowTrace + circshift(swallowTraceBase,[0 round(rand*nT)]);
end
fitpars(3,:) = fitpars(3,:) + swallowMagnitude(1)*swallowTrace;
fitpars(4,:) = fitpars(4,:) + swallowMagnitude(2)*swallowTrace;

% add in random sudden movements in any direction:
suddenTrace = zeros(size(fitpars));
for iS = 1:suddenFrequency
    iT_sudden = ceil(rand*nT);
    suddenTrace(:,iT_sudden:end) = suddenTrace(:,iT_sudden:end) + [suddenMagnitude(1)*((2*rand(3,1))-1); suddenMagnitude(2)*((2*rand(3,1))-1)];
end
fitpars = fitpars+suddenTrace;

%%% uncoment these lines to just have one big rotation
% fitpars = zeros(size(fitpars)); 
% fitpars(6,1:100) = 30;

fitpars = bsxfun(@minus,fitpars,fitpars(:,round(nT/2)));

figure(1)
clf
subplot1(2,1,'Gap',[0 .09],'Max',[.95 1])
s1 = subplot1(1); s2 = subplot1(2);
plotFitPars(fitpars,[],[],[],[s1 s2]);


%% Simulate the effect of that motion

% convert the motion parameters into a set off affine matrices:
fitMats = euler2rmat(fitpars(4:6,:));
fitMats(1:3,4,:) = fitpars(1:3,:); 

% set some things for the recon function:
alignDim = 2; alignIndices = 1:nT; Hxyz = size(rawData); kspaceCentre_xyz = round(Hxyz/2);

% use the recon function just to extract the nufft 'object' st:
[~, st] = applyRetroMC_nufft(rawData,fitMats,alignDim,alignIndices,11,hostVoxDim_mm,Hxyz,kspaceCentre_xyz,-1);
% and use the nufft rather than the nufft_adj function to simulate the rotations:
image_simRotOnly = ifft3s(reshape(nufft(ifft3s(rawData),st),size(rawData)));        
% then apply just the translations:
[~,~,image_simMotion] = applyRetroMC_nufft(fft3s(image_simRotOnly),fitMats,alignDim,alignIndices,11,hostVoxDim_mm,Hxyz,kspaceCentre_xyz,-1);

image_simMotion = ifft3s(image_simMotion);
image_simMotion = image_simMotion / percentile(abs(image_simMotion),95);

% Load both images in a 3D viewer:
SliceBrowser2(cat(4,abs(image_original),abs(image_simMotion)),[0 1.5],{'Original image','Simulated motion'})
set(gcf,'Name','Original image (1) vs Simulated Motion (2)')


%% And how well can this motion be 'undone' again...?

kdata_simMotion = fft3s(image_simMotion);

fitMats_undo = euler2rmat(fitpars(4:6,:)); % keep rotations in the same direction (this will be taken care of in nufft vs nufft_adj)
fitMats_undo(1:3,4,:) = -fitpars(1:3,:);  % swap direction of translations

image_simMoco = applyRetroMC_nufft(kdata_simMotion,fitMats_undo,alignDim,alignIndices,11,hostVoxDim_mm,Hxyz,kspaceCentre_xyz);
image_simMoco = image_simMoco / percentile(abs(image_simMoco),95);

SliceBrowser2(cat(4,abs(image_simMotion),abs(image_simMoco)),[0 1.5],{'Simulated Motion','Simulated Motion-correction'})
set(gcf,'Name','Simulated Motion (1) vs Simulated Motion-correction (2)')


%% Try doing iterative motion-correction (slow on big volumes! - only do this on a sub-sampled dataset or a low-res volume (2mm or less))
 
% tic
% nCGIters = 10; % not sure what the magic number is here... 
% image_simMocoIter = applyRetroMC_nufft(kdata_simMotion,fitMats_undo,alignDim,alignIndices,11,hostVoxDim_mm,Hxyz,kspaceCentre_xyz,nCGIters);
% toc
% image_simMocoIter = image_simMocoIter / percentile(abs(image_simMocoIter),95);
% SliceBrowser2(cat(4,abs(image_simMoco),abs(image_simMocoIter),abs(image_original)),[0 1.5],{'MoCo single NUFFT','MoCo iterative NUFFT','Original Image'})
% set(gcf,'Name','Simulated Motion-correction (1) vs Iterative simulated Motion-correction (2) vs Original Image (3)')
