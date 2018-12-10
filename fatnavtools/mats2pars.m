function pars = mats2pars(Amats)
% Convert affine matrices to rigid-body motion parameters

pars = zeros(6,size(Amats,3));

for iV = 1:size(pars,2)
    thesepars = spm_imatrix(Amats(:,:,iV));
    pars(:,iV) = thesepars(:,1:6);
end
pars(4:6,:) = pars(4:6,:) * 180/pi;