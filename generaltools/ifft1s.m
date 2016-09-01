function out = ifft1s(data, dim, rescale)
% function out = ifft1s(data, dim, rescale)
% 
% Do 1D FFT on dim (default = first) dimension of data and fftshift

if nargin < 3
    rescale = 0;
end

if nargin < 2
    dim = 1;
end

if dim==1 && size(data,1) == 1
    dim = 2;
end

out = fftshift(ifft(ifftshift(data,dim),[],dim),dim);

if rescale
    out = out * sqrt(size(data,dim));
end