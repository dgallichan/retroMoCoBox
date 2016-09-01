function out = get_spm_affmats(P, refVol)

nV = length(P);

if nargin < 2
    refVol = 1;
end

out.pars = zeros(12,nV);
out.mats = zeros(4,4,nV);

for iV = 1:nV
    out.mats(:,:,iV) = P(iV).mat / P(refVol).mat;
    out.pars(:,iV) = spm_imatrix(out.mats(:,:,iV));
end