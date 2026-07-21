function GPUNFT = getGPUNUFFTobjFromCoords(nkx,nky,nkz,dataDims_mm,fullMatrixDims,oversampFactor,sMaps,dcf)



if nargin < 7
    sMaps = [];
end

crds = [nkx(:).'*dataDims_mm(1); nky(:).'*dataDims_mm(2); nkz(:).'*dataDims_mm(3)];

gpu.k = crds/2;

if nargin < 8
    gpu.w = ones(size(gpu.k,2),1);
else
    gpu.w = dcf;
end

gpu.osf = oversampFactor;
gpu.wg = 3;
gpu.sw = 8; % no idea what this should be... 8 is used in example
gpu.imageDim = fullMatrixDims;
gpu.sens = sMaps;

GPUNFT = gpuNUFFT(gpu.k,gpu.w,gpu.osf,gpu.wg,gpu.sw,gpu.imageDim,gpu.sens);


