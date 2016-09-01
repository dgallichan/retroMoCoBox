  function Xk = nufft_table_adj_split(st, X, order, flips, om)
%|function Xk = nufft_table_adj_split(st, X, order, flips, om)
%| adjoint of table-based nufft interpolation.
%| in
%|	st		structure from nufft_init
%|	X [M,nc]	DTFT values (usually nc=1)
%|	order	0|1	0th or 1st-order interpolation
%|			default is 0 for backward compatability
%|	flips	0|1	sign flips? (for real table with even N)
%|	om [M,1]	optional (default st.om)
%| out
%|	Xk [*Kd,nc]	DFT coefficients
%| Copyright 2004-3-30, Jeff Fessler and Yingying Zhang, University of Michigan
%
% 21/4/16 - daniel.gallichan@epfl.ch - try to allow splitting up recon
% across parfor... (3D only)

nSplit = 12; % danielg - no of chunks to split the recon into. I chose 12 because our server has a matlabpool of 12 workers by default.


if nargin < 2, help(mfilename), error(mfilename), end

if ~isvar('order') || isempty(order)
	order = 0; % default 0th order for backward compatability
end

if ~isvar('flips') || isempty(flips)
	flips = zeros(size(st.Nd)); % default no flips for backward compatability
end

if nargin < 4
	om = st.om;
end

dd = length(st.Kd);

if dd ~= 3
    disp('Error, not 3 dimensional - use version without _split.m instead')
    return
end

tm = zeros(size(om));
for id=1:dd
	gam = 2*pi / st.Kd(id);
	tm(:,id) = om(:,id) / gam; % t = omega / gamma
end

if size(X,1) ~= size(om,1)
	error 'X size problem'
end

nc = size(X,2);

% adjoint of phase shift
if isfield(st, 'phase_shift') & ~isempty(st.phase_shift)
	X = X .* repmat(conj(st.phase_shift), [1 nc]);
end

% convert X to complex double for mex file
if ~isa(X, 'double'), X = double(X); end


Xk = zeros(prod(st.Kd), nc,nSplit);

Jd = int32(st.Jd);
Ld = int32(st.Ld);
Kd = int32(st.Kd(1:dd));

Npts = size(X,1);
nPerSplit = ceil(Npts/nSplit);

parfor iSplit = 1:nSplit
    iiSplit = [1:nPerSplit]+(iSplit-1)*nPerSplit;
    iiSplit(iiSplit>Npts) = [];
    arg = {Jd, Ld, tm(iiSplit,:), Kd, int32(order), int32(flips)};
    Xk(:,:,iSplit) = interp3_table_adj_mex(complexify(X(iiSplit,:)), st.h{1}, st.h{2}, st.h{3}, arg{:});
end

Xk = sum(Xk,3);
