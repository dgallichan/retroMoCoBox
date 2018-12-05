function affmats = pars2affmats(pars)
% function affmats = pars2affmats(pars)
% 
% pars is 6 x n matrix of mocopars
% output is 4 x 4 x n affine matrices

npars = size(pars,2);

affmats = zeros(4,4,npars);

for iP = 1:npars
    affmats(:,:,iP) = spm_matrix(pars(:,iP).');
end