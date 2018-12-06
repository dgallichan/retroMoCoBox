function [gdata, st, newData] = applyRetroMC_nufft(kdata_in, alignMats, alignDim, alignIndices, useTable, dataDims_mm, fullMatrixDims, kspaceCentre, cgIterations, oversampFactor)
% The alignMats are a series of 4x4 (rigid-body) matrices describing the
% motion, one for each set of k-space lines with the same motion parameters.
% Displacements (4th column of alignMats) are in mm, or voxels if no 'dataDims_mm' are specified
% 'alignDim' - the k-space dimension along which the motion parameters are available (can be two numbers for unconventional reordering)
% 'alignIndices' - the indices along 'alignDim' for each of the motion matrices. If 'alignDim' is 2D, then this can also be a map of indices
% which were acquired together.
% 'useTable' decides whether to use look-up table for NUFFT gridding - faster but less accurate...
%     - the number becomes 2^n table look-up values
%     - use a negative number to enable parfor in the nufft_table lookup
% 'fullMatrixDims' are the image dimensions after partial Fourier
% 'kspaceCentre' are the centre coordinates of k-space
%
% Requires Jeffrey Fessler's NUFFT code (downloadable from:
% http://web.eecs.umich.edu/~fessler/code/ )
%
% -- Daniel Gallichan, CIBM, EPFL, Lausanne, July 2014
% - 31/7/14 - fixed input default definition problem...
% - 25/6/15 - added option for iterative CG application of the NUFFT
%             (although haven't yet found it to improve any images...!)
% - 17/6/16 - recently I added the option to use 'nufft_table_adj_split.m'
%             which also uses parfor - but this is also a memory hog, so
%             now I've changed it so that only a negative 'useTable' input
%             enables this option
% - 23/8/16 - Added option to change the oversampling factor as the
%             previous default of 1.375 was indeed faster, but perhaps not
%             the best compromise with respect to additional aliasing.
% - 25/8/17 - Added a 'trick' option so that if cgIterations is set to -1
%             then the NUFFT is not applied. Useful if just trying to get
%             at the st object.

if nargin < 10
    oversampFactor = 1.5;
end

if nargin < 9
    cgIterations = 1;
end

if nargin < 8
    kspaceCentre = size(kdata_in)/2 + 1;
end

if nargin < 7
    fullMatrixDims = size(kdata_in);
end

if nargin < 6
    dataDims_mm = size(kdata_in);
end

if nargin < 5
    useTable = 11; 
elseif useTable < 0
    useTable = abs(useTable);
    bUseParfor = 1;
else
    bUseParfor = 0;
end
    

if nargin < 3
    alignDim = 2;
end

if nargin < 4
    alignIndices = 1:size(alignMats,3);
end


if length(alignDim)==1     
    Nt = length(alignIndices);
    maxIndex = max(alignIndices);
    if (maxIndex > size(kdata_in,alignDim))
        disp('Error: indices get too large for the input data along the specified dimension')
        return
    end
else
    % assume that alignIndices is a 2D map along both alignDims...
   idxMap = alignIndices;
   Nt = max(idxMap(:));
   if (Nt~=size(alignMats,3))
       disp('Error: the length of the alignMats should match the indices in idxMap')
       return
   end      
end


Nx = fullMatrixDims(1); Ny = fullMatrixDims(2); Nz = fullMatrixDims(3);
[nx,ny,nz] = size(kdata_in);

[kx, ky, kz] = ndgrid( [1:nx]-kspaceCentre(1), [1:ny]-kspaceCentre(2), [1:nz]-kspaceCentre(3) );
kx = single(2*kx/Nx/dataDims_mm(1));
ky = single(2*ky/Ny/dataDims_mm(2));
kz = single(2*kz/Nz/dataDims_mm(3));


nkx = kx; nky = ky; nkz = kz;

newData = kdata_in;
for iT = 1:Nt
    if length(alignDim)==1     
        switch alignDim
            case 1
                iLine = find(kx==kx(alignIndices(iT),1,1));
            case 2
                iLine = find(ky==ky(1,alignIndices(iT),1));
            case 3
                iLine = find(kz==kz(1,1,alignIndices(iT)));
        end
    else
        switch alignDim(1)
            case 1
                switch alignDim(2)
                    case 2
                        % x + y
                        [iX, iY] = find(idxMap==iT);
                        Npts = length(iX);
                        iX = repmat(iX,[nz 1]);
                        iY = repmat(iY,[nz 1]);
                        iZ = repmat([1:nz],[Npts 1]); iZ = iZ(:);
                        iLine = sub2ind([nx ny nz],iX,iY,iZ);
                    case 3
                        % x + z
                        [iX, iZ] = find(idxMap==iT);
                        Npts = length(iX);
                        iX = repmat(iX,[ny 1]);
                        iZ = repmat(iZ,[ny 1]);
                        iY = repmat([1:ny],[Npts 1]); iY = iY(:);
                        iLine = sub2ind([nx ny nz],iX,iY,iZ);
                end
            case 2
                switch alignDim(2)
                    case 1
                        % y + x
                        [iY, iX] = find(idxMap==iT);
                        Npts = length(iX);
                        iX = repmat(iX,[nz 1]);
                        iY = repmat(iY,[nz 1]);
                        iZ = repmat([1:nz],[Npts 1]); iZ = iZ(:);
                        iLine = sub2ind([nx ny nz],iX,iY,iZ);
                    case 3
                        % y + z
                        [iY, iZ] = find(idxMap==iT);
                        Npts = length(iY);
                        iY = repmat(iY,[nx 1]);
                        iZ = repmat(iZ,[nx 1]);
                        iX = repmat([1:nx],[Npts 1]); iX = iX(:);
                        iLine = sub2ind([nx ny nz],iX,iY,iZ);
                end
            case 3
                switch alignDim(2)
                    case 1
                        % z + x
                        [iZ, iX] = find(idxMap==iT);
                        Npts = length(iX);
                        iX = repmat(iX,[ny 1]);
                        iZ = repmat(iZ,[ny 1]);
                        iY = repmat([1:ny],[Npts 1]); iY = iY(:);
                        iLine = sub2ind([nx ny nz],iX,iY,iZ);
                    case 2
                        % z + y
                        [iY, iZ] = find(idxMap==iT);
                        Npts = length(iY);
                        iY = repmat(iY,[nx 1]);
                        iZ = repmat(iZ,[nx 1]);
                        iX = repmat([1:nx],[Npts 1]); iX = iX(:);
                        iLine = sub2ind([nx ny nz],iX,iY,iZ);
                end
        end
    end    
    newVec = (alignMats(1:3,1:3,iT) * [kx(iLine) ky(iLine) kz(iLine)].' ).';
    nkx(iLine) = newVec(:,1); nky(iLine) = newVec(:,2); nkz(iLine) = newVec(:,3);
    newData(iLine) = newData(iLine).*exp(-1i*pi*( alignMats(1,4,iT)*kx(iLine) + alignMats(2,4,iT)*ky(iLine) + alignMats(3,4,iT)*kz(iLine)));
    if Nt < 1000
        fprintf('.');
    else
        if mod(iT,500)==0
            fprintf('.');
        end
    end
end
fprintf('\n');


crds = [nkx(:).'*dataDims_mm(1); nky(:).'*dataDims_mm(2); nkz(:).'*dataDims_mm(3)];

clear kx ky kz nkx nky nkz

om = [crds.']*pi; % trajectory coords from -pi to pi
Nd = [Nx Ny Nz]; % image dimensions (N1,N2,...,Nd)
Jd = [4 4 4];%	Jd [d]		# of neighbors used (in each direction)
Kd = round(oversampFactor*Nd);%	Kd [d]		FFT sizes (should be >= N1,N2,...)

if useTable
%     st = nufft_init(om, Nd, Jd, Kd,Nd/2,'table',2^useTable,'minmax:kb');
    st = nufft_init(om, Nd, Jd, Kd,Nd/2,'table',2^useTable,'kaiser');
    if bUseParfor
        st.interp_table_adj = @nufft_table_adj_split; % use alternative NUFFT adjoint which uses parfor
    end
else
    st = nufft_init(om, Nd, Jd, Kd,Nd/2,'minmax:kb');
end

clear crds kdata_in om

if cgIterations==1    
    gdata = nufft_adj_single(newData(:),st);
else
    if cgIterations==-1
        % trick to skip applying the nufft if we just want to get at the st
        % object
        gdata = [];
    else
        disp('...Using CG recon version....')
        E = newOperator(@(x) reshape(nufft(reshape(x,Nd),st),[],1), @(x) reshape(nufft_adj_single(reshape(x,[],1),st),[],1));
        gdata = reshape(CGbasic(E,newData(:),'maxIters',cgIterations),Nd);
    end
end

 