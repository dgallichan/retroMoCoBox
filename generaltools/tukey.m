function W = tukey(n,kc,w)
% function W = tukey(n,kc,w)
%
% Tukey windowing function, as described on p496 of Bernstein2004

if nargin < 3, w = 0.15; end
if nargin < 2, kc = 0.85; end

if length(n)==1
    kvals = linspace(-1,1,n);
else
    kvals = n;
end

W = cos((pi*(abs(kvals)-kc))/(2*w)).^2;
W(abs(kvals) < kc) = 1;
W(abs(kvals) > (kc+w)) = 0;



