function newAffMats = recentre_affmats(affMats,iCentreVol)

nV = size(affMats,3);
newAffMats = zeros(size(affMats));

for iV = 1:nV
    newAffMats(:,:,iV) = affMats(:,:,iV) / affMats(:,:,iCentreVol);    
end