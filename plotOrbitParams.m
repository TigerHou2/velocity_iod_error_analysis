%% plotOrbitParams.m
%
% Author:
%   Tiger Hou
%
% Description:
%   Generates plots for velocity initial orbit determination (VIOD) error
%   with respect to variations in orbital parameters, including semi-major
%   axis, eccentricity, and true anomaly.
%
% Notes:
%   In the case of varying orbital parameters, we choose to hold spacecraft
%   parameters constant. This means all simulations will be conducted with
%   equal noise, number of measurements, and duration of observation.
%   Specifically, we will use the following:
%       - Noise = 1e-6 DU/TU
%       - Number of measurements = 3
%       - Duration of observation = 0.1 orbital period

%% initialization

close all hidden
clear;clc;init

% load orbital parameter config case studies
[smaVect,eccVect,taVect,noise,names] = ...
    load_orbit_cases();

% define gravitational parameter as 1 in canonical units
mu = 1; % DU^3/TU^2

% define orbital parameters
SMA  = 1e6; % DU
ECC  = 0.1; % nd
INC  = 24;  % deg
RAAN = 30;  % deg
AOP  = 30;  % deg
TA   = 0;   % deg

% convert to radians and combine parameters
orbitParams = [SMA,ECC,deg2rad([INC,RAAN,AOP,TA])];

% set the index of the measurement at which error will be compared
compIdx = 1;

% spacecraft parameters (fixed)
numObsv = 20; % nd, number of measurements
duration = 0.1; % nd, measurement duration as fraction of the orbit period

% Monte Carlo settings
numSims = 3000;
rngSeed = 1;

% sanity checks
if compIdx > numObsv
    error(['The comparison index must be less or equal to '...
           'the number of measurements!'])
end

% data logging
filepath = 'data\orbitParams.mat';
meanOrig = nan(length(eccVect),length(taVect),length(smaVect));
meanHodo = nan(length(eccVect),length(taVect),length(smaVect));
 stdOrig = nan(length(eccVect),length(taVect),length(smaVect));
 stdHodo = nan(length(eccVect),length(taVect),length(smaVect));

%% generate ground truth data and perturbations

% initialize cell array of velocity measurements and position references
vArray = cell(length(smaVect),length(eccVect),length(taVect));
rArray = cell(length(smaVect),length(eccVect),length(taVect));
% populate velocity and position cell arrays
for i = 1:length(smaVect)
    SMA = smaVect(i);
for j = 1:length(eccVect)
    ECC = eccVect(j);
for k = 1:length(taVect)
    TA = taVect(k);
    % load new case
    orbitParams(1) = SMA;
    orbitParams(2) = ECC;
    orbitParams(6) = TA;
    % initialize empty velocity matrix
    v = nan(numObsv,3);
    % calculate mean anomaly
    period = 2*pi*sqrt(SMA^3/mu);
    E = 2 * atan(sqrt((1-ECC)/(1+ECC))*tan(TA/2));
    M = E - ECC*sin(E);
    
    % determine measurement locations
    %------ Option 1: uniform TRUE anomaly distribution ------
%         Me = M + duration*2*pi;
%         Ee = kepler(Me,ECC);
%         fe = 2 * atan(sqrt((1+ECC)/(1-ECC))*tan(Ee/2));
%         if fe < TA
%             fe = fe + 2*pi;
%         end
%         fVect = linspace(TA,fe,numObsv);
    %------ Option 2: uniform MEAN anomaly distribution ------
        Mvect = M + linspace(0,duration,numObsv) * 2*pi;
        Evect = kepler(Mvect,ECC);
        fVect = 2 * atan(sqrt((1+ECC)/(1-ECC))*tan(Evect/2));

    % fetch all measurements
    for p = 1:numObsv
        orbitParams(6) = fVect(p);
        [~,v(p,:)] = Get_Orb_Vects(orbitParams,mu);
        if p == compIdx
            rRef = Get_Orb_Vects(orbitParams,mu);
        end
    end
    vArray{i,j,k} = v;
    rArray{i,j,k} = rRef;
end %taVect
end %eccVect
end %smaVect

% initialize cell array of noise values
nArray = cell(1,numSims);
% populate noise cell array
rng(rngSeed);
for n = 1:numSims
    noiseVect = randn(numObsv,3);
    noiseVect = noiseVect ./ vecnorm(noiseVect,2,2) * noise;
    nArray{n} = noiseVect;
end %numSims

%% apply perturbations and perform VIOD

% line formats
origWidth = 1.5;
hodoWidth = 0.3;
origFormat = '-x';
hodoFormat = '--k';

% color formats
cmap = hsv(length(smaVect));

% initialize cell array for custom labeling
labelArray = cell(length(smaVect),length(eccVect),2);

for i = 1:length(smaVect)
    SMA = smaVect(i);
for j = 1:length(eccVect)
    ECC = eccVect(j);
    
% initialize position error matrix
errOrig = nan(numSims,length(taVect));
errHodo = nan(numSims,length(taVect));

% perform Monte Carlo sim
for n = 1:numSims
for k = 1:length(taVect)
    v = vArray{i,j,k};
    noiseVect = nArray{n};
    rRef = rArray{i,j,k};
    % do orbit determination
    % --- note the scaling dor the original method: this is because
    % --- precision issues arise when using canonical units. 
%     rOrig = viod((v+noiseVect)*1e4,mu*1e12)/1e4;
    rOrig = hodo(v+noiseVect,mu);
    rHodo = hodoHyp(v+noiseVect,mu);
    % we choose to compare the position estimate at the first measurement
    rOrig = rOrig(compIdx,:)';
    rHodo = rHodo(compIdx,:)';
    % store error data
    errOrig(n,k) = norm(rOrig-rRef) / norm(rRef) * 100;
    errHodo(n,k) = norm(rHodo-rRef) / norm(rRef) * 100;
end %taVect
end %numSims

% plot results for each case
color = cmap(i,:);
taDeg = rad2deg(taVect);

% --- mean error
figure(1)
latexify(30,18)
hold on
plot(taDeg,mean(errOrig),origFormat,'Color',color,'LineWidth',origWidth)
plot(taDeg,mean(errHodo),hodoFormat,'LineWidth',hodoWidth)
hold off
legend('Hodograph Method','Hyperfit Method','Location','SouthEast')
xlabel('True Anomaly, deg')
ylabel('Position Error Avg. \%')
% store annotation data point
labelArray{i,j,1} = [taDeg(end),mean(errOrig(:,end))];
latexify(18)

% --- standard deviation
figure(2)
latexify(30,18)
hold on
plot(taDeg,std(errOrig),origFormat,'Color',color,'LineWidth',origWidth)
plot(taDeg,std(errHodo),hodoFormat,'LineWidth',hodoWidth)
hold off
legend('Hodograph Method','Hyperfit Method','Location','SouthEast')
xlabel('True Anomaly, deg')
ylabel('Position Error StDev \%')
% store annotation data point
labelArray{i,j,2} = [taDeg(end),std(errOrig(:,end))];
latexify(18)

% --- log data
meanOrig(j,:,i) = mean(errOrig)' / 100;
meanHodo(j,:,i) = mean(errHodo)' / 100;
 stdOrig(j,:,i) =  std(errOrig)' / 100;
 stdHodo(j,:,i) =  std(errHodo)' / 100;

end %eccVect
end %smaVect

% save data to file
save(filepath,'meanOrig','meanHodo','stdOrig','stdHodo',...
              'smaVect','eccVect','taVect');

% print data in latex table formatting
regexprep(...
regexprep(...
regexprep(latex(vpa(sym(meanHodo(:,:,1)),10)),...
    '([0-9]+\.[0-9]+)','${num2str(str2num($1),''%.3e'')}'),...
    '(e[+-][0-9]+)',' \\times 10^{${strtok($1,''e'')}}'),...
    '((?<=\{[+-])[0-9])','${regexprep($1,''^0*'','''')}')
regexprep(...
regexprep(...
regexprep(latex(vpa(sym(meanHodo(:,:,2)),10)),...
    '([0-9]+\.[0-9]+)','${num2str(str2num($1),''%.3e'')}'),...
    '(e[+-][0-9]+)',' \\times 10^{${strtok($1,''e'')}}'),...
    '((?<=\{[+-])[0-9])','${regexprep($1,''^0*'','''')}')
regexprep(...
regexprep(...
regexprep(latex(vpa(sym(meanHodo(:,:,3)),10)),...
    '([0-9]+\.[0-9]+)','${num2str(str2num($1),''%.3e'')}'),...
    '(e[+-][0-9]+)',' \\times 10^{${strtok($1,''e'')}}'),...
    '((?<=\{[+-])[0-9])','${regexprep($1,''^0*'','''')}')

%% add annotations

figure(1)
set(gca, 'YScale', 'log')
setgrid
expand(0.05,0.10,0.02,0.03)
axis tight
axlim = axis(gca) .* [1 1 0.3 3];
xlim(axlim(1:2))
ylim(axlim(3:4))
for i = 1:length(smaVect)
    SMA = smaVect(i);
for j = 1:length(eccVect)
    ECC = eccVect(j);
    p = labelArray{i,j,1};
    % convert from log scale to linear scale
    p(2) = (axlim(4)-axlim(3)) * log(p(2)/axlim(3)) / log(axlim(4)/axlim(3));
    [figx,figy] = ax2fig(p(1),p(2));
    annotation('textbox',[figx figy-0.08 1 .1],...
               'String',['e = ' num2str(ECC,4)],...
               'EdgeColor','none',...
               'Interpreter','latex',...
               'FontSize',14)
end
end

figure(2)
set(gca, 'YScale', 'log')
setgrid
expand(0.05,0.10,0.02,0.03)
axis tight
axlim = axis(gca) .* [1 1 0.3 3];
xlim(axlim(1:2))
ylim(axlim(3:4))
for i = 1:length(smaVect)
    SMA = smaVect(i);
for j = 1:length(eccVect)
    ECC = eccVect(j);
    p = labelArray{i,j,2};
    % convert from log scale to linear scale
    p(2) = (axlim(4)-axlim(3)) * log(p(2)/axlim(3)) / log(axlim(4)/axlim(3));
    [figx,figy] = ax2fig(p(1),p(2));
    annotation('textbox',[figx figy-0.08 1 .1],...
               'String',['e = ' num2str(ECC,4)],...
               'EdgeColor','none',...
               'Interpreter','latex',...
               'FontSize',14)
end
end

%% add orbit SMA legend

figure(1)
ax = gca;
ax2 = copyobj(ax,gcf);
delete(get(ax2,'Children')) % delete its children
hold(ax2,'on')
for i = 1:length(smaVect)
    h(i) = plot(sin(1:10),'Parent',ax2,...
                          'Color',cmap(i,:),...
                          'LineWidth',1,...
                          'DisplayName',['SMA = ' num2str(smaVect(i),4)]);
end
hold(ax2,'off')
legend(ax2, 'Location','NorthWest')
for i = 1:length(smaVect)
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
for i = 1:length(smaVect)
    h(i) = plot(sin(1:10),'Parent',ax2,...
                          'Color',cmap(i,:),...
                          'LineWidth',1,...
                          'DisplayName',['SMA = ' num2str(smaVect(i),4)]);
end
hold(ax2,'off')
legend(ax2, 'Location','NorthWest')
for i = 1:length(smaVect)
    set(h(i),'XData',[],'YData',[]);
end
set(ax2,'Color','none',...
    'XTick',[],'YTick',[],...
    'Box','off') % make legend axes transparent
ax2.XLabel.String = '';
ax2.YLabel.String = '';