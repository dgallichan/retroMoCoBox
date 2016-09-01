function out = GrappaReco3D_arb(sig,kernel,ws,stepSize,startPos)
% function out = GrappaReco3D_arb(sig,kernel,ws,stepSize,startPos)
%
% Daniel Gallichan, December 2012. Small modification to code from
% bernd.jung@uniklinik-freiburg.de to allow arbitrary GRAPPA kernel
% and multiple 'startPos' (i.e. the corner of the first kernel which will
% be moved by 'stepSize' to reconstruct the data).
%
% e.g if your data look like this:
%
% . X . X .
% X . X . X 
% . X . X . 
% X . X . X
% . X . X . 
%
%
% you would want to use a stepSize of [2 2], but you would need two
% starting positions to fill in all the missing data. Specify this with
% startPos = [1 1; 2 2]; (each row is considered a pair of starting
% positions)
%
%
% Note that I also removed the capability of including points in the third
% dimension (readout for 3D acceleration) because it keeps the code simpler
% and in some (albeit not very extensive...) tests it didn't actually make
% any difference.
%


% Determination of kData-size
[ny,nz,nc]=size(sig);

kernelsize_y = size(kernel,1);
kernelsize_z = size(kernel,2);
[ySrcVc,zSrcVc]= find(kernel==1);
[yTargVc,zTargVc]= find(kernel==0.5);

noOfSrcPoints = size(ySrcVc,1);
noOfTargPoints = size(yTargVc,1);

% create a border with zeros around sig-matrix
yBorder = ceil(kernelsize_y/stepSize(1))*stepSize(1); zBorder = ceil(kernelsize_z/stepSize(2))*stepSize(2);
% yBorder = 0; zBorder=0;
nyExt = ny+2*yBorder;
nzExt = nz+2*zBorder;
sigIn = zeros( nyExt, nzExt, nc);
sigIn(yBorder+1:end-yBorder, 1+zBorder:end-zBorder, :) = sig;
	
% do cyclic boundary conditions (but still excludes corners!)
if yBorder>0
    sigIn(1:yBorder,1+zBorder:end-zBorder,:) = sig(end-yBorder+1:end,:,:);
    sigIn(end-yBorder+1:end,1+zBorder:end-zBorder,:) = sig(1:yBorder,:,:);
end
if zBorder>0
    sigIn(1+yBorder:end-yBorder,1:zBorder,:) = sig(:,end-zBorder+1:end,:);
    sigIn(1+yBorder:end-yBorder,end-zBorder+1:end,:) = sig(:,1:zBorder,:);
end

sigOut = zeros(size(sigIn));

    
for iStart = 1:size(startPos,1)     
    for zind= startPos(iStart,2):stepSize(2):nzExt-kernelsize_z+1
        for yind = startPos(iStart,1):stepSize(1):nyExt-kernelsize_y+1
            tmpMx = sigIn(yind:yind+kernelsize_y-1,zind:zind+kernelsize_z-1,:);
            tmp_srcMx = complex(zeros(noOfSrcPoints,size(tmpMx,3)));
            for i=1:noOfSrcPoints
                tmp_srcMx(i,:) = tmpMx(ySrcVc(i),zSrcVc(i),:);
            end
            src = reshape(tmp_srcMx, 1, nc*noOfSrcPoints);
            tmpMx = reshape(src*ws,[noOfTargPoints nc]);
            for i=1:noOfTargPoints
                sigOut(yind-1+yTargVc(i),zind-1+zTargVc(i),:) = tmpMx(i,:);
            end                      
        end
    end
end

% put the original data back in place in the output
out = sig;

% now add in the new values where they have been calculated (overwriting if
% they overlap)
sigOut = sigOut( yBorder+1:yBorder+ny, zBorder+1:zBorder+nz, :);
mask_sigOut = zeros(size(sigOut));
mask_sigOut(sigOut~=0) = 1;
out = out.*(1-mask_sigOut) + sigOut;
