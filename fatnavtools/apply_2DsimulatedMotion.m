function [kdata_simMotion, st] = apply_2DsimulatedMotion(kdata_in,mpars)

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
oversampFactor = 2; % keep this high (i.e. around 2) for forward simulation of artifacts

alignMats_swapTrans = alignMats;
alignMats_swapTrans(1:3,4,:) = -alignMats_swapTrans(1:3,4,:);


% run once just to be able to extract the 'st' object
[~, st] = applyRetroMC_nufft_2D(kdata_in, alignMats_swapTrans, alignDim, alignIndices, useTable, dataDims_mm, fullMatrixDims, kspaceCentre,-1,oversampFactor);

% and use the nufft rather than the nufft_adj function to simulate the rotations:
kdata_simRotOnly = reshape(nufft(ifft2s(kdata_in),st),size(kdata_in));        
% then apply just the translations:
[~,~,kdata_simMotion] = applyRetroMC_nufft_2D(kdata_simRotOnly,alignMats_swapTrans,alignDim,alignIndices,useTable,dataDims_mm,fullMatrixDims,kspaceCentre,-1,oversampFactor);
