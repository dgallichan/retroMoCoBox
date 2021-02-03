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

    
    %%% The following code was in place for many years, but I'm not sure
    %%% where I got it from originally...! It seems to work correctly until
    %%% double oblique slices are used, then it deviates from expected
    %%% behaviour (non-normalised rotation matrix...)
%     out(:,:,iT) = [  ct*cps cp*sps+sp*st*cps sp*sps-cp*st*cps;
%         -ct*sps cp*cps-sp*st*sps sp*cps-cp*st*sps;
%         st        -sp*ct           cp*ct;];

%%% instead use the method from spm_matrix: 
% (note that sign on st in R2 is swapped compared to spm_matrix.m so that
%  euler2rmat.m and rotmat2euler.m remain a compatible pair)

R1 = [1    0    0;
      0   cp   sp; 
      0  -sp   cp];
  
R2 = [ct  0   -st;
      0    1   0  ;
      st 0   ct];
  
R3 = [cps sps 0;
     -sps cps 0;
      0   0 1];
  
  out(:,:,iT) = R1*R2*R3;

    
end