function [gdata, st, newData] = applyRetroMC_nufft_2D(kdata_in, alignMats, alignDim, alignIndices, useTable, dataDims_mm, fullMatrixDims, kspaceCentre)
% rehash of main function to work on simulated 2D data
%
% This is a bit 'quick-and-dirty' implementation, but as 2D is so much
% faster than 3D, it's probably fine for the purposes of exploring
% motion-corruption and quality metrics, etc.

% run([getmatlabpath() '/../matlabdownloads/fessler/setup.m']);

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


Nx = fullMatrixDims(1); Ny = fullMatrixDims(2);
[nx,ny] = size(kdata_in);

[kx, ky] = ndgrid( [1:nx]-kspaceCentre(1), [1:ny]-kspaceCentre(2) );
kx = single(2*kx/Nx/dataDims_mm(1));
ky = single(2*ky/Ny/dataDims_mm(2));


nkx = kx; nky = ky; 

newData = kdata_in;
for iT = 1:Nt
    if length(alignDim)==1     
        switch alignDim
            case 1
                iLine = find(kx==kx(alignIndices(iT),1,1));
            case 2
                iLine = find(ky==ky(1,alignIndices(iT),1));
        end
    else
       disp('Error, 2D alignDim not supported for 2D code...')
       return
    end    
    newVec = (alignMats(1:2,1:2,iT) * [kx(iLine) ky(iLine)].' ).';
    nkx(iLine) = newVec(:,1); nky(iLine) = newVec(:,2); 
    newData(iLine) = newData(iLine).*exp(-1i*pi*( alignMats(1,4,iT)*kx(iLine) + alignMats(2,4,iT)*ky(iLine)));
    fprintf('.');
end
fprintf('\n');


crds = [nkx(:).'*dataDims_mm(1); nky(:).'*dataDims_mm(2)];

clear kx ky nkx nky

om = [crds.']*pi; % trajectory coords from -pi to pi
Nd = [Nx Ny]; % image dimensions (N1,N2,...,Nd)
Jd = [4 4];%	Jd [d]		# of neighbors used (in each direction)
Kd = 2*Nd;%	Kd [d]		FFT sizes (should be >= N1,N2,...)

if useTable
    st = nufft_init(om, Nd, Jd, Kd,Nd/2,'table',2^useTable,'minmax:kb');
else
    st = nufft_init(om, Nd, Jd, Kd,Nd/2,'minmax:kb');
end

clear crds kdata_in om

gdata = nufft_adj_single(newData(:),st);

 