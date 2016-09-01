function out = fft1s(data, dim, rescale)
% function out = fft1s(data, dim, rescale)
% 
% Do 1D FFT on dim (default=1) dimension of data and fftshift

if nargin < 3
    rescale = 0;
end

if nargin < 2
    dim = 1;
end

if dim==1 && size(data,1) == 1
    dim = 2;
end

out = fftshift(fft(ifftshift(data,dim),[],dim),dim);

if rescale
    out = out / sqrt(size(data,dim));
end

