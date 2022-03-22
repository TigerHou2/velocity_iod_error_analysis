%% plotSpcraftParams.m
%
% Author:
%   Tiger Hou
%
% Description:
%   Generates plots for velocity initial orbit determination (VIOD) error
%   with respect to variations in spacecraft parameters, including
%   measurement noise and measurement scheduling.
%
% Notes:
%   In the case of varying spacecraft parameters, we choose to hold orbital
%   parameters constant. This means all simulations will be conducted with
%   equal semi-major axis, eccentricity, and true anomaly. The orbital
%   parameters chosen are as follows:
%       - SMA = 1e5 DU
%       - ECC = 0.5
%       - TA  = 20 deg

%% initialization

close all hidden
clear;clc;init

% load spacecraft parameter config case studies
[noiseVect,durVect,obsVect] = load_spcraft_cases();

% define gravitational parameter as 1 in canonical units
mu = 1; % DU^3/TU^2

% define orbital parameters
SMA  = 1.56e5; % DU
ECC  = 0.5; % nd
INC  = 24;  % deg
RAAN = 30;  % deg
AOP  = 30;  % deg
TA   = 0;  % deg

% convert to radians and combine parameters
orbitParams = [SMA,ECC,deg2rad([INC,RAAN,AOP,TA])];

% Monte Carlo settings
numSims = 3000;
rngSeed = 1;

%% iterate through all cases

% line formats
origWidth = 1.5;
hodoWidth = 0.3;
origFormat = '-x';
hodoFormat = '--k';

% color formats
cmap = hsv(length(durVect));

% initialize cell array for custom labeling
labelArray = cell(length(durVect),length(obsVect),2);

for i = 1:length(durVect)
    duration = durVect(i);
for j = 1:length(obsVect)
    numObsv = obsVect(j);
    
% initialize position error matrix
errOrig = nan(numSims,length(noiseVect));
errHodo = nan(numSims,length(noiseVect));
% reset random number generator
rng(rngSeed);
% calculate mean anomaly
period = 2*pi*sqrt(SMA^3/mu);
e = orbitParams(2);
f = orbitParams(6);
E = 2 * atan(sqrt((1-e)/(1+e))*tan(f/2));
M = E - e*sin(E);
% determine measurement locations
Mvect = M + linspace(0,duration,numObsv) * 2*pi;
Evect = kepler(Mvect,e);
fVect = 2 * atan(sqrt((1+e)/(1-e))*tan(Evect/2));
% store ground truth position vector at starting true anomaly
rRef = Get_Orb_Vects(orbitParams,mu);

% perform Monte Carlo sim
for n = 1:numSims
noiseBase = randn(numObsv,3);
noiseBase = noiseBase ./ vecnorm(noiseBase,2,2) ;

for k = 1:length(noiseVect)
    % load new case
    noise = noiseBase * noiseVect(k);
    % initialize empty velocity matrix
    v = nan(numObsv,3);
    % fetch all measurements
    for p = 1:numObsv
        orbitParams(6) = fVect(p);
        [~,v(p,:)] = Get_Orb_Vects(orbitParams,mu);
    end
    % restore true anomaly in orbit parameter
    orbitParams(6) = f;
    % do orbit determination
    % --- note the scaling for the original method: this is because
    % --- precision issues arise when using canonical units. 
%     rOrig = viod((v+noise)*1e4,mu*1e12)/1e4;
    rOrig = hodo(v+noise,mu);
    rHodo = hodoHyp(v+noise,mu);
    % we choose to compare the position estimate at the first measurement
    rOrig = rOrig(1,:)';
    rHodo = rHodo(1,:)';
    % store error data
    errOrig(n,k) = norm(rOrig-rRef) / norm(rRef);
    errHodo(n,k) = norm(rHodo-rRef) / norm(rRef);
end %nVect

end %numSims

% plot results for each case
color = cmap(i,:);

% plot results for each case
% --- mean error
figure(1)
latexify(24,18)
noiseScale = 1e-6;
hold on
plot(noiseVect/noiseScale,...
        mean(errOrig),origFormat,'LineWidth',origWidth,'Color',color)
plot(noiseVect/noiseScale,...
        mean(errHodo),hodoFormat,'LineWidth',hodoWidth)
hold off
legend('Energy Method','Hodograph Method','Location','SouthEast')
xlabel(['Noise, ' num2str(noiseScale) ' DU/TU'])
ylabel('Position Error Avg., fraction of SMA')
latexify(18)
% store annotation data point
labelArray{i,j,1} = [noiseVect(end)/noiseScale,mean(errOrig(:,end))/SMA];

% --- standard deviation
figure(2)
latexify(24,18)
noiseScale = 1e-6;
hold on
plot(noiseVect/noiseScale,...
        std(errOrig),origFormat,'LineWidth',origWidth,'Color',color)
plot(noiseVect/noiseScale,...
        std(errHodo),hodoFormat,'LineWidth',hodoWidth)
hold off
legend('Energy Method','Hodograph Method','Location','SouthEast')
xlabel(['Noise, ' num2str(noiseScale) ' DU/TU'])
ylabel('Position Error StDev, fraction of SMA')
latexify(18)
% store annotation data point
labelArray{i,j,2} = [noiseVect(end)/noiseScale,std(errOrig(:,end))/SMA];

end %obsVect
end %durVect

%% add annotations

figure(1)
set(gca, 'YScale', 'log')
setgrid
expand(0.07,0.12,0.07,0.05)
axis tight
axlim = axis(gca) .* [1 1 0.3 3];
xlim(axlim(1:2))
ylim(axlim(3:4))
for i = 1:length(durVect)
    duration = durVect(i);
for j = 1:length(obsVect)
    numObsv = obsVect(j);
    p = labelArray{i,j,1};
    % convert from log scale to linear scale
    p(2) = (axlim(4)-axlim(3)) * log(p(2)/axlim(3)) / log(axlim(4)/axlim(3));
    [figx,figy] = ax2fig(p(1),p(2));
    annotation('textbox',[figx figy-0.08 1 .1],...
               'String',['obsv = ' num2str(numObsv)],...
               'EdgeColor','none',...
               'Interpreter','latex',...
               'FontSize',14)
end
end

figure(2)
set(gca, 'YScale', 'log')
setgrid
expand(0.07,0.12,0.07,0.05)
axis tight
axlim = axis(gca) .* [1 1 0.3 3];
xlim(axlim(1:2))
ylim(axlim(3:4))
for i = 1:length(durVect)
    duration = durVect(i);
for j = 1:length(obsVect)
    numObsv = obsVect(j);
    p = labelArray{i,j,2};
    % convert from log scale to linear scale
    p(2) = (axlim(4)-axlim(3)) * log(p(2)/axlim(3)) / log(axlim(4)/axlim(3));
    [figx,figy] = ax2fig(p(1),p(2));
    annotation('textbox',[figx figy-0.08 1 .1],...
               'String',['obsv = ' num2str(numObsv)],...
               'EdgeColor','none',...
               'Interpreter','latex',...
               'FontSize',14)
end
end

%% add measurement duration legend

figure(1)
ax = gca;
ax2 = copyobj(ax,gcf);
delete(get(ax2,'Children')) % delete its children
hold(ax2,'on')
for i = 1:length(durVect)
    h(i) = plot(sin(1:10),'Parent',ax2,...
                          'Color',cmap(i,:),...
                          'LineWidth',1,...
                          'DisplayName',['Dur = ' num2str(durVect(i),4)]);
end
hold(ax2,'off')
legend(ax2, 'Location','NorthWest')
for i = 1:length(durVect)
    set(h(i),'XData',[],'YData',[]);
end
set(ax2,'Color','none',...
    'XTick',[],'YTick',[],...
    'Box','off') % make legend axes transparent
ax2.XLabel.String = '';
ax2.YLabel.String = '';

figure(2)
ax = gca;
ax2 = copyobj(ax,gcf);
delete(get(ax2,'Children')) % delete its children
hold(ax2,'on')
for i = 1:length(durVect)
    h(i) = plot(sin(1:10),'Parent',ax2,...
                          'Color',cmap(i,:),...
                          'LineWidth',1,...
                          'DisplayName',['Dur = ' num2str(durVect(i),4)]);
end
hold(ax2,'off')
legend(ax2, 'Location','NorthWest')
for i = 1:length(durVect)
    set(h(i),'XData',[],'YData',[]);
end
set(ax2,'Color','none',...
    'XTick',[],'YTick',[],...
    'Box','off') % make legend axes transparent
ax2.XLabel.String = '';
ax2.YLabel.String = '';