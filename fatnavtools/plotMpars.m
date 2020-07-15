function plotMpars(Amats,t,iRef)
% Plot the effect of a series of affine matrices on 3 points located at
% (x/y/z) = 50mm

nT = size(Amats,3);

if nargin < 2
    t = 1:nT;
end

if nargin < 3
%     iRef = round(nT/2);
    iRef = 1;
end

pts = [1 0 0; 0 1 0; 0 0 1]*50;
pts(:,4) = 1;

Vpts = zeros(nT,3,3);
for iT = 1:nT
    for iV = 1:3        
        thisV = Amats(:,:,iT)*pts(iV,:).';
        Vpts(iT,iV,:) = thisV(1:3);
    end
end

Vdisp = bsxfun(@minus,Vpts,Vpts(iRef,:,:));
% Vdisp2 = squeeze(ssos(Vdisp,2)); % this looks at movement of each point in all 3 dimensions
Vdisp2 = zeros(nT,3);
Vdisp2(:,1) = squeeze(ssos(Vdisp(:,[2 3],1),2));
Vdisp2(:,2) = squeeze(ssos(Vdisp(:,[1 3],2),2));
Vdisp2(:,3) = squeeze(ssos(Vdisp(:,[1 2],3),2));

figure
set(gcf,'Position',[        1962         670        1083         318])
plot(t,Vdisp2)
xlim([t(1) t(end)])
xlabel('Time')
ylabel('Displacements (mm)')
grid on 
grid minor
fontScale(1.4)
