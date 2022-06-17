% run_retroMoCoDemo_FatNavRecon.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Load in some example FatNav data from a real experiment (2mm resolution,
% 4x4 GRAPPA acceleration) and do the reconstruction 
%
% Can be used to get a feeling for the process involved in reconstructing
% each FatNav if you don't have any full datasets of your own to hand.
% 
% Example data can be downloaded from here: https://bit.ly/retroData_4 (186 MB) 
% 
% Data consists of the ACS prescan (auto-calibration data, fully sampled
% centre of k-space) along with 15 volumes of 2mm data acquired with 4x4
% undersampling in both PE directions.

run('addRetroMoCoBoxToPath.m')

load ../exampleData/exampleFatNavdata_15fatnavs.mat

%%

FatNavRes_mm = 2;

Ax = 4; Ay = 4; % Hard-coded for a 'known' navigator...
nc_orig = 32;

% apply shift due to fat offset frequency (this assumes BW of 1950 Hz, and a fat/water shift of 3.5ppm)
offsetsXYZ = [0 0 298*3.5/1950]; 

nT = 15;

dataSize = size(ACSdata);

nx = 176/FatNavRes_mm;
ny = 256/FatNavRes_mm;
nz = 256/FatNavRes_mm;
nc = size(ACSdata,2);

ACSdata = squeeze(ACSdata);         
ACSdata = permute(ACSdata,[4 3 1 2]);
ACSdata(nx,ny,nz,1) = 0; % extend to full size
ACSdata = circshift(ACSdata,double(round([nx/2-dataSize(4)/2 ny/2-dataSize(3)/2 0 0])));  
% no idea why it should be necessary to put 'double' here - but I got the error 'invalid shift type: must be a finite, nonsparse, real integer vector

[ox, oy, oz] = ndgrid(linspace(-1,1,nx),linspace(-1,1,ny),linspace(-1,1,nz));
offsetPhase = exp(1i*pi*( offsetsXYZ(1)*ox + offsetsXYZ(2)*oy + offsetsXYZ(3)*oz) );
ACSdata = bsxfun(@times,ACSdata,offsetPhase);

ACSdata = ifft1s(ACSdata,3);
ACSdata = ACSdata(:,:,end:-1:1,:);

ACSims = ifft2s(ACSdata);
ACSim = ssos(ACSims);

ovACS = orthoview(ACSim);
set(gcf,'Name','Low-res image from the ACS calibration scan')

ACSdata = zeroCrop(ACSdata); % after making an image at full resolution, drop back to just where the data really is


   
%% calculate the GRAPPA weights separately for all 'virtual slices' of the 3D volume

gx = 2; gy = 2; % size of GRAPPA kernel
grapKernel = zeros(gx + (gx-1)*(Ax-1), gy + (gy-1) * (Ay-1));

% source points are marked with 1
grapKernel(1:Ax:end, 1:Ay:end) = 1;

% target points are marked with 0.5
if gx == 1,    startx = 2; else  startx = 2 + (floor(gx/2) - 1) * Ax; end
if gy == 1,    starty = 2; else  starty = 2 + (floor(gy/2) - 1) * Ay; end

grapKernel(startx:startx+Ax-2, starty:starty+Ay-2) = 0.5;
grapKernel(startx:startx+Ax-1, starty:starty+Ay-2) = 0.5;
grapKernel(startx:startx+Ax-2, starty:starty+Ay-1) = 0.5;

stepSize = [Ax Ay];
all_grapW = zeros(gx*gy*nc,((Ax*Ay)-1)*nc,128);
zMax = nz; % don't need to go outside the brain... would save time to detect the top of the brain here - but could be difficult to make robust...

tic
disp('... Calculating GRAPPA weights ....')
for iZ = 1:zMax
    
    thisACS = squeeze(ACSdata(:,:,iZ,:));
    all_grapW(:,:,iZ) = GrappaCalib3D_arb(thisACS,grapKernel,0,1);
    
    fprintf('.');
end
fprintf('\n');
disp('...Done...')
timingReport.calculateGRAPPAweights = toc;
    
    
%% Do reconstruction of each individual FatNav
    
nxMeas = 61;
nyMeas = 93;
startPos = [nx-nxMeas+1 ny-nyMeas+1];

timingEachFatNav  = zeros(nT,1);
timeBeforeFatNavs = clock;

fatnav_volumes = zeros(88,128,128,15);

for iT = 1 % can go up to 15 with this example data, but might take a while...
    
    full_GRAPPArecons = int16(zeros(nx,ny,nz));
    
%     fatnavdata = twix_obj.FatNav(:,:,:,:,1,1,1,1,iT);
    fatnavdata = zeros(128,32,nyMeas,nxMeas);    
    fatnavdata(:,:,1:4:end,1:4:end) = fatnavdata_15(:,:,:,:,iT);
    
    fatnavdata = permute(fatnavdata,[4 3 1 2]);
    fatnavdata(nx,ny,nz,1) = 0; % extend to full size
    
    fatnavdata = circshift(fatnavdata,double([nx-nxMeas ny-nyMeas 0 0]));
    
    if any(offsetsXYZ)
        fatnavdata = bsxfun(@times,fatnavdata,offsetPhase);
    end
    
    fatnavdata = ifft1s(fatnavdata,3);
    fatnavdata = fatnavdata(:,:,end:-1:1,:);
    
    tic
    for iZ = 1:zMax
        thisfatnavdata = squeeze(fatnavdata(:,:,iZ,:));
        thisfatnavdata(nx,ny,1) = 0;
        
        grapK = GrappaReco3D_arb(thisfatnavdata,grapKernel,all_grapW(:,:,iZ),stepSize,startPos);
        
        % apply a Tukey window k-space filter - it's not really clear what
        % the 'optimal' parameters for this would be, but as the images are
        % going to be used for co-registration, a small amount of smoothing
        % seems to be appropriate.
        f1 = tukey(nx,.7,.3); f2 = tukey(ny,.8,.3); 
        grapK = bsxfun(@times,grapK,f1.');
        grapK = bsxfun(@times,grapK,f2);
        
        grapIms = ifft2s(grapK);
        grapIm = ssos(grapIms)*nx*ny*nz;
        full_GRAPPArecons(:,:,iZ) = int16(4095*grapIm/20); % use 20 as first signal to clip...
        
    end
    
    fatnav_volumes(:,:,:,iT) = full_GRAPPArecons;
        
      
    ov1 = orthoview(ssos(ifft2s(fatnavdata)),'drawIms',0);
    ov2 = orthoview(full_GRAPPArecons,'drawIms',0);
   
    figure(1)
    clf
    set(gcf,'Position',[ 500   172   858   926]);    
    subplot1(2,3); 
    subplot1(1); imab(ov1.im1); subplot1(2); imab(ov1.im2); subplot1(3); imab(ov1.im3);
    subplot1(4); imab(ov2.im1); subplot1(5); imab(ov2.im2); subplot1(6); imab(ov2.im3);
    subplot1(2)
    title('Direct FFT of 4x4 accelerated 3D-FatNav')
    subplot1(5)
    title('GRAPPA reconstruction of the same data')
    subplot1(1); ylabel(['FatNav ' num2str(iT) ' out of ' num2str(nT)]);
    subplot1(4); ylabel(['FatNav ' num2str(iT) ' out of ' num2str(nT)]);    
    drawnow
        
end

    





