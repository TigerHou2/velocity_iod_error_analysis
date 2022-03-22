%% speed.m
%
% Author:
%   Tiger Hou
%
% Description:
%   Compares the speed of three VIOD methods.
%

%% initialization
close all
clear;clc
addpath('..\fcns_od')
addpath('..\fcns_orb')
addpath('..\fcns_vis')
addpath('..\fcns_misc')
latexify

%% setup

% orbit definition
mu = 1;
a = 1e5;
e = 0.5;
i = deg2rad(0);
o = deg2rad(0);
w = deg2rad(0);
f = deg2rad(160);
orbitParams = [a,e,i,o,w,f];

% observations
dM = 0.1*2*pi;
obsVect = ceil(linspace(3,55,8));

% simulation
numSims = 1000;
selObsv = 1;
timeDat = nan(3,length(obsVect));

% noise
noise = 3e-6;

% line styles
EN = 'rs:'; % energy method
HD = 'bo--'; % hodograph method
IM = 'kd-.';  % improved method
format = {EN,HD,IM};

%% vary number of measurements

for i = 1:length(obsVect)
    numObsv = obsVect(i);
    
    % noise
    nGauss = normrnd(0,noise,numObsv,1,numSims);
    nGauss = repmat(nGauss,1,3,1);
    ncube = randn(numObsv,3,numSims);
    ncube = ncube ./ vecnorm(ncube,2,2);
    ncube = ncube .* nGauss;
    
    errDat = nan(numSims,length(obsVect),3);
    v = nan(numObsv,3);

    % find measurement positions by true anomaly
    E0 = 2 * atan(sqrt((1-e)/(1+e))*tan(f/2));
    M0 = E0 - e*sin(E0);
    M = M0 + dM;
    Mvect = linspace(M0,M,numObsv);
    Evect = kepler(Mvect,e);
    fvect = 2 * atan(sqrt((1+e)/(1-e))*tan(Evect/2));
    fvect = mod(fvect,2*pi);
    fend = fvect(end);

    % ground truth data
    for j = 1:numObsv
        orbitParams(6) = fvect(j);
        [~,v(j,:)] = Get_Orb_Vects(orbitParams,mu);
    end

    % get position reference
    orbitParams(6) = fvect(selObsv);
    rRef = Get_Orb_Vects(orbitParams,mu);
    
    % Monte Carlo
    tic
    for s = 1:numSims
        n = ncube(:,:,s);
        r_EN = viod((v+n)*1e4,mu*1e12)/1e4;
        r_EN = r_EN(selObsv,:)';
        r_HD = hodo(v+n,mu);
        errDat(s,i,1) = norm(r_EN-rRef) / norm(rRef) * 100;
    end
    timeDat(1,i) = toc;
    tic
    for s = 1:numSims
        n = ncube(:,:,s);
        r_HD = hodo(v+n,mu);
        r_HD = r_HD(selObsv,:)';
        errDat(s,i,2) = norm(r_HD-rRef) / norm(rRef) * 100;
    end
    timeDat(2,i) = toc;
    tic
    for s = 1:numSims
        n = ncube(:,:,s);
        r_IM = hodoHyp(v+n,mu);
        r_IM = r_IM(selObsv,:)';
        errDat(s,i,3) = norm(r_IM-rRef) / norm(rRef) * 100;
    end
    timeDat(3,i) = toc;
end

figure(1)
hold on
for i = 1:3
    xVar = obsVect;
    yVar = timeDat(i,:);
    plot(xVar,yVar,format{i},'MarkerFaceColor','none','LineWidth',1.25);
end
hold off
set(gca,'YScale','log')
xlabel('Number of Measurements')
ylabel('Computation Time, s')
legend('Energy','Hodograph','Improved','Location','Best')
latexify(15,15,24)
setgrid
expand
print('speed','-depsc')

%% tabulate data
tab = latex(vpa(sym([obsVect',timeDat']),3)); % 3 sig-figs
tab = regexprep(tab,'\\\',' \\\\ \\hline ${newline}'); % outline + return
tab = regexprep(tab,'c\}','c\} ${newline}'); % return for first line
tab = regexprep(tab,'\\end\{','\\\\ \\hline ${newline} \\end\{'); % return for last line
tab = regexprep(tab,'([0-9]+\.[0-9]+)','${strcat(($1),'' s'')}'); % add units
tab = regexprep(tab,'.0 s',''); % format first column to obsv #
disp(tab)