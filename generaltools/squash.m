function [xr, mask] = squash(x, m)

% function [ret, mask] = squash(x)
%
% performs ret = reshape(x,prod(size(x)),1);
%
% [ret, mask] = squash(x, m)
% masks and then squashes, returns the mask
% m can be a mask or a threshold
%
% see unsquash
   
[xr, dims]=sq(x);
   
if(nargin > 1)
if(size(m)==1),
   coords = find(sq(xr>m));
   xr = xr(coords);
   mask = thresh(x,m);
else
   coords = find(sq(m)>0);
   xr = xr(coords);
   mask = m;
end;
else
    mask = ones(dims);
end;
   
function [ret, dims] = sq(x)
   
dims = size(x);
ret = reshape(x,prod(size(x)),1);