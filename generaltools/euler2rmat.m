function out = euler2rmat(phi, theta, psi)

if nargin==1
    eMat = phi;
    nT = size(eMat,2);
else
    nT = 1;
end
    
out = zeros(3,3,nT);

for iT = 1:nT
    
    if nT > 1
        phi = eMat(1,iT)*pi/180; theta = eMat(2,iT)*pi/180; psi = eMat(3,iT)*pi/180;
    else
        phi = phi*pi/180; theta = theta*pi/180;  psi = psi*pi/180;
    end
    
    sp = sin(phi);    cp = cos(phi);
    st = sin(theta);  ct = cos(theta);
    sps = sin(psi);   cps = cos(psi);
    
    out(:,:,iT) = [  ct*cps cp*sps+sp*st*cps sp*sps-cp*st*cps;
        -ct*sps cp*cps-sp*st*sps sp*cps-cp*st*sps;
        st        -sp*ct           cp*ct;];
    
end