function [spMat_src, spMat_targ] = GrappaCalib3D_arb_idx(acsDims,kernel)
%
%
% Daniel Gallichan: modification of Bernd Jung's code for GRAPPA - to
% specify manually the kernel
%
% - 5/2/14 - changed the definition of lambda so that it is relative to
% norm of the ACS data
%
% - 11/3/16 - change default when no regularization is set to use pinv
%             instead of left-divide, which seems to be more sensible
%             default behaviour.
%             Also added option to not actually calculate the weights, but
%             only the src and targ matrices so that they can be
%             concatenated if desired in external functions.
%
% get just the indices for a particular kernel and acs dimensions



% Determination of acs kData-size
% [nyacs,nzacs,nc]=size(acs);
nyacs = acsDims(1); nzacs = acsDims(2); nc = acsDims(3);

% bx2=fix(bx/2); % remove readout direction capability for now...

kernelsize_y = size(kernel,1);
kernelsize_z = size(kernel,2);
[ySrcVc,zSrcVc]= find(kernel==1);
[yTargVc,zTargVc]= find(kernel==0.5);

noOfSrcPoints = size(ySrcVc,1);
noOfTargPoints = size(yTargVc,1);

% src  = complex(zeros( (nzacs-kernelsize_z)*(nyacs-kernelsize_y),nc*noOfSrcPoints ));
% targ = complex(zeros( (nzacs-kernelsize_z)*(nyacs-kernelsize_y),nc*noOfTargPoints ));

ssz = (nzacs-kernelsize_z+1)*(nyacs-kernelsize_y+1);
spMat_src = logical(sparse( ssz*noOfSrcPoints, nzacs*nyacs ));
spMat_targ = logical(sparse( ssz*noOfTargPoints, nzacs*nyacs ));

idxMat = reshape(1:nzacs*nyacs,[nyacs nzacs]);

% Fitting procedure
cnt=0;
for zind=1:nzacs-kernelsize_z+1
    for yind=1:nyacs-kernelsize_y+1
%         for xind= bx2+1 : nxacs-bx2
            cnt=cnt+1;
%             tmpMx = acs(yind:yind+kernelsize_y-1,zind:zind+kernelsize_z-1,:);
            tmpMx_idx = idxMat(yind:yind+kernelsize_y-1,zind:zind+kernelsize_z-1);
%             tmp_srcMx = complex(zeros(noOfSrcPoints,size(tmpMx,3)));
            tmp_srcMx_idx = zeros(noOfSrcPoints,1);
            for i=1:noOfSrcPoints
%                 tmp_srcMx(i,:) = tmpMx(ySrcVc(i),zSrcVc(i),:);
                tmp_srcMx_idx(i) = tmpMx_idx(ySrcVc(i),zSrcVc(i));
                spMat_src(cnt+(i-1)*ssz,tmp_srcMx_idx(i)) = true;
            end
%             src(cnt,:)  = reshape(tmp_srcMx, 1, nc*noOfSrcPoints);
                                       
                        
%             tmp_targMx = complex(zeros(noOfTargPoints,size(tmpMx,3)));
            tmp_targMx_idx = zeros(noOfTargPoints,1);
            for i=1:noOfTargPoints
%                 tmp_targMx(i,:) = tmpMx(yTargVc(i),zTargVc(i),:);
                tmp_targMx_idx(i) = tmpMx_idx(yTargVc(i),zTargVc(i));
                spMat_targ(cnt+(i-1)*ssz,tmp_targMx_idx(i)) = true;
            end
%             targ(cnt,:) = reshape(tmp_targMx, 1, nc*noOfTargPoints);
            
%             tmp_targMx_idx = tmp_targMx_idx + (0:noOfTargPoints-1)'*ssz;
            
%         end
    end
end

% 
% if (bDoInversion)
%     
%     % danielg - add tikhonov regularisation, but not sure if it correct for
%     % multicolumn 'targ'... changes result, but not necessarily in a good
%     % way...
%     if (lambda)
%         lambda = lambda*norm(src(:));
%         %     ws = [src; lambda*eye(size(src,2))] \ [targ; zeros(size(src,2),size(targ,2))];
%         ws = pinv(src'*src + lambda^2*eye(size(src,2))) * src' * targ;
%     else
%         %     ws = src \ targ;
%         ws = pinv(src)*targ;
%         
%     end
%     
% else
%     ws = [];
% end




