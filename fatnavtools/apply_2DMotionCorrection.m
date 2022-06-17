function [im_corrected, st, newData] = apply_2DMotionCorrection(kdata_in,mpars,cgIters,oversampFactor)

if nargin < 4
    oversampFactor = 2;
end

if nargin < 3
    cgIters = 1;
end

Nres = size(mpars,2);
alignMats = zeros(4,4,Nres);
this_mpars = mpars;
for iT = 1:Nres
    alignMats(1:3,1:3,iT) = euler2rmat(0,0,this_mpars(3,iT));
    alignMats(1:2,4,iT) = this_mpars(1:2,iT);
end
alignMats(4,4,:) = 1;
alignDim = 2;
alignIndices = 1:Nres;
useTable = 0;
dataDims_mm = [1 1 1];
fullMatrixDims = [Nres Nres 1];
kspaceCentre = [Nres/2+1 Nres/2+1 1];

% alignMats_swapTrans = alignMats;
% alignMats_swapTrans(1:3,4,:) = -alignMats_swapTrans(1:3,4,:); % reversal
% only for forward simulation of artifacts

[im_corrected, st, newData] = applyRetroMC_nufft_2D(kdata_in, alignMats, alignDim, alignIndices, useTable, dataDims_mm, fullMatrixDims, kspaceCentre,cgIters,oversampFactor);

   
