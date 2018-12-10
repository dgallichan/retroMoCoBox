function out = orthoview(vol,varargin)
% function orthoview(vol,...)

[cvox voxdim clims useColorbar useMip useMeanip useTitles drawIms interpFact useNewFig] = process_options(varargin,...
    'centre',[],'voxdim',[],'clims',[],'colorbar',0,'mip',0,'meanip',0,'useTitles',1,'drawIms',1,'interpFact',1,'useNewFig',1);

if ~any(isreal(vol(:)))
    vol = abs(vol);
    disp('Data are complex, only using magnitude...')
end

if size(vol,4) > 1
    disp('Data are 4 dimensional - only showing first volume...')
    vol = vol(:,:,:,1);
end

if isempty(cvox)
   cvox = round(size(vol)/2);
end

if ~isempty(voxdim)
  vdims = abs(voxdim);
  vdims = vdims + (vdims==0)*1;  % get rid of zero dimensions
else
  vdims = [1 1 1];
end

if ~isempty(clims)
  irange = [clims(1) clims(2)];
else
  % Use the robust min and max
  irange = percentile(vol,[2 98]);
end

if ( abs(irange(2) - irange(1))==0),
  % Use absolute range if robust is no good
  irange = [min(min3(vol)) max(max3(vol))];
end
if ( abs(irange(2) - irange(1))==0),
  irange = [irange(1) irange(1)+1];
end

out.vdims = vdims;
out.irange = irange;

if useMip && useMeanip
    disp('Error - MIP and MeanIP are mutually exclusive...!')
    return
end

if ~useMip && ~useMeanip
    out.im1 = vol(:,:,cvox(3));
    out.im2 = squeeze(vol(:,cvox(2),:));
    out.im3 = squeeze(vol(cvox(1),:,:));
elseif useMip
    out.im1 = max(vol,[],3);
    out.im2 = squeeze(max(vol,[],2));
    out.im3 = squeeze(max(vol,[],1));
else % now useMeanip must be true!
    out.im1 = mean(vol,3);
    out.im2 = squeeze(mean(vol,2));
    out.im3 = squeeze(mean(vol,1));
end

if drawIms
    if useNewFig
        figure
        set(gcf,'Position',[    50   720   950  340])
    else
        clf
    end
    subplot1(1,3)
    subplot1(1)
    imab(gs(out.im1,0,interpFact),irange(:)')
    set(gca,'DataAspectRatio',[vdims(2) vdims(1) 1]);
    if useTitles, title('xy'); end;
    subplot1(2)
    imab(gs(out.im2,0,interpFact),irange(:)')
    set(gca,'DataAspectRatio',[vdims(3) vdims(1) 1]);
    if useTitles, title('xz'); end;
    subplot1(3)
    imab(gs(out.im3,0,interpFact),irange(:)')
    set(gca,'DataAspectRatio',[vdims(3) vdims(2) 1]);
    if useTitles, title('yz'); end;
    colormap(gray)
    if useColorbar
        hc = colorbar;
        set(hc,'Position',[    0.9389    0.159    0.0200    0.726])
    end
end

matDims = size(vol);
maxDim = max(matDims);
ims = NaN(maxDim,maxDim,3);
ims(1:matDims(1),1:matDims(2),1) = out.im1;
ims(1:matDims(1),1:matDims(3),2) = out.im2;
ims(1:matDims(2),1:matDims(3),3) = out.im3;
shifts = round((matDims-maxDim)/2);
ims(:,:,1) = circshift(ims(:,:,1),-[shifts(1) shifts(2)]);
ims(:,:,2) = circshift(ims(:,:,2),-[shifts(1) shifts(3)]);
ims(:,:,3) = circshift(ims(:,:,3),-[shifts(2) shifts(3)]);
out.ims = ims;

x1 = size(out.im1,1)+size(out.im2,1)+size(out.im3,1);
y1 = max([size(out.im1,2) size(out.im2,2) size(out.im3,2)]);
out.oneIm = zeros(x1,y1);
out.oneIm(1:size(out.im1,1),1:size(out.im1,2)) = out.im1;
out.oneIm(size(out.im1,1)+[1:size(out.im2,1)],1:size(out.im2,2)) = out.im2;
out.oneIm(size(out.im1,1)+size(out.im2,1)+[1:size(out.im3,1)],1:size(out.im3,2)) = out.im3;
