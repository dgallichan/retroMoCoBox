function out = fft2s(data,useRescale)
% function out = fft2s(data,useRescale)
% 
% Do 2D FFT on each slice of data and fftshift

if nargin < 2
    useRescale = 0;
end

out = fftshift(fftshift(fft2(ifftshift(ifftshift(data,1),2)),1),2); 
if useRescale
    nx = size(out,1); ny = size(out,2);
    out = out / sqrt(nx) / sqrt(ny);
end