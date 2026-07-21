function pars = mats2pars(Amats, doInv)
% Convert affine matrices to rigid-body motion parameters

if nargin < 2
    doInv = 0;
end
    

pars = zeros(6,size(Amats,3));

for iV = 1:size(pars,2)
    if doInv==0
        thesepars = spm_imatrix(Amats(:,:,iV));
    elseif doInv==1
        thesepars = spm_imatrix(inv(Amats(:,:,iV)));
    elseif doInv==2
        thisMat = Amats(:,:,iV);
        thisMat(1:3,1:3) = thisMat(1:3,1:3).';
        thisMat(1:3,4) = -thisMat(1:3,4);
        thesepars = spm_imatrix(thisMat);
    end
    pars(:,iV) = thesepars(:,1:6);
end
pars(4:6,:) = pars(4:6,:) * 180/pi;