function ha = plotFitPars(fitPars,t,climsXYZ,climsRTP,isVertical)
% function ha = plotFitPars(fitPars,t,climsXYZ,climsRTP,isVertical)
%
% 16/3/16 - Added a little trick so that 'isVertical' can also be of the
% form [handle1 handle2] for the two subplots to plot into

if size(fitPars,1)==4 && size(fitPars,2)==4 % assume affine matrices
    Amats = fitPars;
    % match the way the get_spm_affmats function was written
    fitPars = mats2pars(Amats);   
end

if nargin < 5
    isVertical = 1;
end

if nargin < 4 || isempty(climsRTP)
    maxRot = max(abs(squash(fitPars(4:6,:))));
    if maxRot > 0
        climsRTP = 1.5*maxRot*[-1 1];
    else
        climsRTP = [-1 1];
    end    
end
if length(climsRTP)==1
    climsRTP = climsRTP*[-1 1];
end

if nargin < 3 || isempty(climsXYZ)
    maxDisp = max(abs(squash(fitPars(1:3,:))));
    if maxDisp > 0
        climsXYZ = 1.5*maxDisp*[-1 1];
    else
        climsXYZ = [-1 1];
    end
end
if length(climsXYZ)==1
    climsXYZ = climsXYZ*[-1 1];
end

if nargin < 2 || isempty(t)
    t = 1:size(fitPars,2);
    timeUnits = '';
else
    if isduration(t)
        timeUnits = '';
        t.Format = 'hh:mm:ss';
    else
    timeUnits = '(mins)';
end
end



lstyle = {'linewidth',2};

%%

cOrder =    [      0         0    1.0000;...
         0    0.5000         0;...
    1.0000         0         0;...
         0    0.7500    0.7500;...
    0.7500         0    0.7500;...
    0.7500    0.7500         0;...
    0.2500    0.2500    0.2500];

if isnumeric(isVertical)
    figure
end

displacements = sqrt(sum(fitPars(1:3,:).^2,1));
RMS_displacement = sqrt(mean(displacements.^2));

rotations = sqrt(sum(fitPars(4:6,:).^2,1));
RMS_rot = sqrt(mean(rotations.^2));

if isVertical(1)==1
    set(gcf,'Pos',[              1810         498         726         480])
    subplot1(2,1,'Gap',[0 .09],'Max',[.95 1])
    ha.s1 = subplot1(1);
    ha.s2 = subplot1(2);
elseif isVertical(1)==0
    %%% Horizontal:
    set(gcf,'Pos',[        1810         637        1534         341])
    subplot1(1,2,'YTickL','All','Gap',[.05 0],'Min',[0.05 .2],'Max', [1 .9])
    ha.s1 = subplot1(1);
    ha.s2 = subplot1(2);
else
    ha.s1 = isVertical(1);
    ha.s2 = isVertical(2);
end

subplot(ha.s1)
set(gca,'ColorOrder',cOrder);
plot(t,fitPars(1:3,:).',lstyle{:})
% axis([0 t(end) climsXYZ])
set(gca,'XLim',[t(1) t(end)]);
set(gca,'YLim',climsXYZ);
grid on
ylabel('Displacements (mm)')
legend('x','y','z','location','north','orientation','horizontal','autoupdate','off')
if isVertical(1)==0 || ~isnumeric(isVertical)
    xlabel(['Time ' timeUnits])
end
title(['RMS displacement: ' num2str(RMS_displacement,'%.2f') ' mm'])
subplot(ha.s2)
set(gca,'ColorOrder',cOrder);
plot(t,fitPars(4:6,:).',lstyle{:})
% axis([0 t(end) climsRTP])
set(gca,'XLim',[t(1) t(end)]);
set(gca,'YLim',climsRTP);
title(['RMS rotations: ' num2str(RMS_rot,'%.2f') ' degrees'])
grid on
xlabel(['Time ' timeUnits])
ylabel('Rotations (degrees)')
% legend('\phi','\theta','\psi','location','northwest')
legend('Pitch','Roll','Yaw','location','north','orientation','horizontal','autoupdate','off')

if verLessThan('matlab','8.4')
    fontScale(1.6)
end


