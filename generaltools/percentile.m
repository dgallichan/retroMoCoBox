function [y] = percentile(x,p)
% y = percentile(x,p)
%   p is a number between 0 and 100 (the percentile value)
%   x is the data vector (it will reshape any matrix into a single column)
%
% I know there is the function prctile.m which does this, but it seems
% stupid to check out a license for a toolbox just for this simple
% function!

% if (p>100), p=100; end
% if (p<0),   p=0;   end
p(p>100)=100; p(p<0) = 0;
n = numel(x);
xx = sort(x(:));
y = xx(1 + round((n-1)*p/100))';
