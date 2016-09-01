function ha = plotFitPars(fitPars,t,climsXYZ,climsRTP,isVertical)
% function ha = plotFitPars(fitPars,t,climsXYZ,climsRTP,isVertical)
%
% 16/3/16 - Added a little trick so that 'isVertical' can also be of the
% form [handle1 handle2] for the two subplots to plot into

if nargin < 5
    isVertical = 1;
end

if nargin < 4 || isempty(climsRTP)
    climsRTP = 1.5*max(abs(squash(fitPars(4:6,:))))*[-1 1];
end
if length(climsRTP)==1
    climsRTP = climsRTP*[-1 1];
end

if nargin < 3 || isempty(climsXYZ)
    climsXYZ = 1.5*max(abs(squash(fitPars(1:3,:))))*[-1 1];
end
if length(climsXYZ)==1
    climsXYZ = climsXYZ*[-1 1];
end

if nargin < 2 || isempty(t)
    t = 1:size(fitPars,2);
    timeUnits = '';
else
    timeUnits = '(mins)';
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
    s1 = subplot1(1);
    s2 = subplot1(2);
elseif isVertical(1)==0
    %%% Horizontal:
    set(gcf,'Pos',[        1810         637        1534         341])
    subplot1(1,2,'YTickL','All','Gap',[.05 0],'Min',[0.05 .2],'Max', [1 .9])
    s1 = subplot1(1);
    s2 = subplot1(2);
else
    s1 = isVertical(1);
    s2 = isVertical(2);
end

subplot(s1)
set(gca,'ColorOrder',cOrder);
plot(t,fitPars(1:3,:).',lstyle{:})
axis([0 t(end) climsXYZ])
grid on
ylabel('Displacements (mm)')
legend('x','y','z','location','north','orientation','horizontal')
if isVertical(1)==0 || ~isnumeric(isVertical)
    xlabel(['Time ' timeUnits])
end
title(['RMS displacement: ' num2str(RMS_displacement,'%.2f') ' mm'])
subplot(s2)
set(gca,'ColorOrder',cOrder);
plot(t,fitPars(4:6,:).',lstyle{:})
axis([0 t(end) climsRTP])
title(['RMS rotations: ' num2str(RMS_rot,'%.2f') ' degrees'])
grid on
xlabel(['Time ' timeUnits])
ylabel('Rotations (degrees)')
% legend('\phi','\theta','\psi','location','northwest')
legend('Pitch','Roll','Yaw','location','north','orientation','horizontal')

if verLessThan('matlab','8.4')
    fontScale(1.6)
end


