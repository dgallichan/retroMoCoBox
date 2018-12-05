function outMats = moveFrame(inMats,applyMat)

nT = size(inMats,3);
outMats = zeros(4,4,nT);

for iT = 1:nT
    outMats(:,:,iT) = applyMat*inMats(:,:,iT)/applyMat;
end