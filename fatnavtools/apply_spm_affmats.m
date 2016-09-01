function Vout = apply_spm_affmats(V,Mats)

nV = length(V);

if (nV)~=size(Mats,3)
    disp('')
    disp('Error: there should be a 4x4 matrix for the n volumes from V')
    disp('')
    return
end

Vout = V;

for iV = 1:nV
    Vout(iV).mat = Mats(:,:,iV)*Vout(1).mat;
end