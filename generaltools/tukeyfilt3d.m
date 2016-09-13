function out = tukeyfilt3d(data,kc,w,isRadial)
% function W = tukey(n,kc,w,isRadial)
%
% Tukey windowing function, as described on p496 of Bernstein2004

if nargin < 4, isRadial = 0; end
if nargin < 3, w = 0.15; end
if nargin < 2, kc = 0.85; end


nx = size(data,1); ny = size(data,2); nz = size(data,3);
sz = size(data);

if ~isRadial
    Wx = tukey(nx,kc,w)';
    Wy = tukey(ny,kc,w);
    Wz = reshape(tukey(nz,kc,w),[1 1 nz]);
    Wxyz = bsxfun(@times,bsxfun(@times,Wx,Wy),Wz);
else
    [kx, ky, kz] = ndgrid(linspace(-1,1,nx),linspace(-1,1,ny),linspace(-1,1,nz));
    kvals = sqrt(kx.^2+ky.^2+kz.^2);
    Wxyz = cos((pi*(abs(kvals)-kc))/(2*w)).^2;
    Wxyz(abs(kvals) < kc) = 1;
    Wxyz(abs(kvals) > (kc+w)) = 0;
end

Wsz = sz; Wsz(1:3) = 1;
out = repmat(Wxyz,Wsz).*data;