function st = getNUFFTobjFromCoords(nkx,nky,nkz,dataDims_mm,fullMatrixDims,oversampFactor,useTable,bUseParfor)

if nargin < 8
    bUseParfor = 0;
end

crds = [nkx(:).'*dataDims_mm(1); nky(:).'*dataDims_mm(2); nkz(:).'*dataDims_mm(3)];

om = [crds.']*pi; % trajectory coords from -pi to pi
% Nd = [Nx Ny Nz]; % image dimensions (N1,N2,...,Nd)
Nd = fullMatrixDims;
Jd = [4 4 4];%	Jd [d]		# of neighbors used (in each direction)
Kd = round(oversampFactor*Nd);%	Kd [d]		FFT sizes (should be >= N1,N2,...)

if useTable
    st = nufft_init(om, Nd, Jd, Kd,Nd/2,'table',2^useTable,'kaiser');
    if bUseParfor
        st.interp_table_adj = @nufft_table_adj_split; % use alternative NUFFT adjoint which uses parfor
    end
else
    st = nufft_init(om, Nd, Jd, Kd,Nd/2,'minmax:kb');
end


