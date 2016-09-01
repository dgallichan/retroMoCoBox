function out = ifft3s(data, useRescale)
% function out = ifft3s(data, useRescale)
% 
% Do 3D FFT on first 3 dims of data and fftshift


if nargin < 2
    useRescale = 0;
end

out = ifft1s(ifft2s(data),3);

if useRescale
    nx = size(out,1); ny = size(out,2); nz = size(out,3);
    out = out * sqrt(nx) * sqrt(ny) * sqrt(nz);
end