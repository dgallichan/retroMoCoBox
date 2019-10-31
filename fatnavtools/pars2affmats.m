function affmats = pars2affmats(pars)
% function affmats = pars2affmats(pars)
% 
% pars is 6 x n matrix of mocopars (3 translations, 3 rotations in degrees)
% output is 4 x 4 x n affine matrices

npars = size(pars,2);

pars(4:6,:) = pars(4:6,:) * pi/180; % make into radians for SPM

affmats = zeros(4,4,npars);

for iP = 1:npars
    affmats(:,:,iP) = spm_matrix(pars(:,iP).');
end