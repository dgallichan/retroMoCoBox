function FD = calc_framewiseDisplacements(Amats)
% Calculation 'framewise displacement' motion scores, based on 
% Tisdall et al, MRM 2012

nT = size(Amats,3);

FD = zeros(nT-1,1);

for iT = 2:nT
    
    % get difference in affine matrices:
    thisDiff = Amats(:,:,iT) / Amats(:,:,iT-1);
    mpars = mats2pars(thisDiff); % dx dy dz (mm) r1 r2 r3 (degrees)
    
    dx = mpars(1); dy = mpars(2); dz = mpars(3);
    cx = cosd(mpars(4)); cy = cosd(mpars(5)); cz = cosd(mpars(6));
    sx = sind(mpars(4)); sy = sind(mpars(5)); sz = sind(mpars(6));
    
    % get net angle of rotation (degrees)    
    thisTheta = abs(acosd(.5*(-1 + cx*cy + cx*cz + cy*cz + sx*sy*sz)));
    
    % get 'DeltaR' (Tisdall Eq2)
    DeltaR = 64 * sqrt( (1 - cosd(thisTheta))^2 + sind(thisTheta)^2 );
    
    FD(iT-1) = DeltaR + sqrt(dx^2+dy^2+dz^2);
    
end