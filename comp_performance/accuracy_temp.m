%% accuracy.m
%
% Author:
%   Tiger Hou
%
% Description:
%   Compares the accuracy of three VIOD methods.
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
f = deg2rad(90);
orbitParams = [a,e,i,o,w,f];

% observations
numObsv = 20;
dM = 0.5*2*pi;

% simulation
numSims = 3000;
selObsv = 1;

% noise
noise = 3e-6;
ncube = randn(numObsv,3,numSims);
ncube = ncube ./ vecnorm(ncube,2,2);

% line styles
EN = 'rs:'; % energy method
HD = 'bo--'; % hodograph method
IM = 'kd-.';  % improved method
format = {EN,HD,IM};

%% vary noise

nVect = linspace(0,noise*4,9);

errDat = nan(numSims,length(nVect),3);
v = nan(numObsv,3);

for i = 1:length(nVect)
    noise_loc = nVect(i);
    % prepare measurement noise
    nGauss = normrnd(0,noise_loc,numObsv,1,numSims);
    nGauss = repmat(nGauss,1,3,1);
    ncube_loc = ncube .* nGauss;
    
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
    for s = 1:numSims
        n = ncube_loc(:,:,s);
        r_EN = hodoSuper(v+n,mu);
        r_EN = r_EN(selObsv,:)';
        r_HD = hodo(v+n,mu);
        r_HD = r_HD(selObsv,:)';
        r_IM = hodoHyp(v+n,mu);
        r_IM = r_IM(selObsv,:)';
        errDat(s,i,1) = norm(r_EN-rRef) / norm(rRef) * 100;
        errDat(s,i,2) = norm(r_HD-rRef) / norm(rRef) * 100;
        errDat(s,i,3) = norm(r_IM-rRef) / norm(rRef) * 100;
    end
end

figure(1)
hold on
for i = 1:3
    xVar = nVect;
    yVar = mean(errDat(:,:,i));
    yStd =  std(errDat(:,:,i));
    plot(xVar,yVar,format{i},'MarkerFaceColor','none','LineWidth',1.25);
end
hold off
xlabel('1-$\sigma$ Noise, $\frac{DU}{TU}$')
ylabel('RMSE$(\tilde{\mathbf{r}})$, \%')
legend('Super','Hodograph','Improved','Location','NorthWest')
latexify(12,12,22)
setgrid
expand
% print('accuracy_noise','-dpdf','-bestfit')

% restore orbit
orbitParams = [a,e,i,o,w,f];

%% vary eccentricity

eccVect = linspace(0,0.9,10);

errDat = nan(numSims,length(eccVect),3);
v = nan(numObsv,3);

% prepare measurement noise
nGauss = normrnd(0,noise,numObsv,1,numSims);
nGauss = repmat(nGauss,1,3,1);
ncube_loc = ncube .* nGauss;

for i = 1:length(eccVect)
    e_loc = eccVect(i);
    orbitParams(2) = e_loc;
    
    % find measurement positions by true anomaly
    E0 = 2 * atan(sqrt((1-e_loc)/(1+e_loc))*tan(f/2));
    M0 = E0 - e_loc*sin(E0);
    M = M0 + dM;
    Mvect = linspace(M0,M,numObsv);
    Evect = kepler(Mvect,e_loc);
    fvect = 2 * atan(sqrt((1+e_loc)/(1-e_loc))*tan(Evect/2));
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
    for s = 1:numSims
        n = ncube_loc(:,:,s);
        r_EN = hodoSuper(v+n,mu);
        r_EN = r_EN(selObsv,:)';
        r_HD = hodo(v+n,mu);
        r_HD = r_HD(selObsv,:)';
        r_IM = hodoHyp(v+n,mu);
        r_IM = r_IM(selObsv,:)';
        errDat(s,i,1) = norm(r_EN-rRef) / norm(rRef) * 100;
        errDat(s,i,2) = norm(r_HD-rRef) / norm(rRef) * 100;
        errDat(s,i,3) = norm(r_IM-rRef) / norm(rRef) * 100;
    end
end

figure(2)
hold on
for i = 1:3
    xVar = eccVect;
    yVar = mean(errDat(:,:,i));
    yStd =  std(errDat(:,:,i));
    plot(xVar,yVar,format{i},'MarkerFaceColor','none','LineWidth',1.25);
end
hold off
xlabel('Eccentricity')
ylabel('RMSE$(\tilde{\mathbf{r}})$, \%')
% legend('Energy','Hodograph','Improved','Location','NorthWest')
latexify(12,12,22)
setgrid
expand
% print('accuracy_ecc','-dpdf','-bestfit')

% restore orbit
orbitParams = [a,e,i,o,w,f];

%% vary true anomaly

taVect = deg2rad(linspace(0,330,12));

errDat = nan(numSims,length(taVect),3);
v = nan(numObsv,3);

% prepare measurement noise
nGauss = normrnd(0,noise,numObsv,1,numSims);
nGauss = repmat(nGauss,1,3,1);
ncube_loc = ncube .* nGauss;

for i = 1:length(taVect)
    f_loc = taVect(i);
    orbitParams(6) = f_loc;
    
    % find measurement positions by true anomaly
    E0 = 2 * atan(sqrt((1-e)/(1+e))*tan(f_loc/2));
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
    for s = 1:numSims
        n = ncube_loc(:,:,s);
        r_EN = hodoSuper(v+n,mu);
        r_EN = r_EN(selObsv,:)';
        r_HD = hodo(v+n,mu);
        r_HD = r_HD(selObsv,:)';
        r_IM = hodoHyp(v+n,mu);
        r_IM = r_IM(selObsv,:)';
        errDat(s,i,1) = norm(r_EN-rRef) / norm(rRef) * 100;
        errDat(s,i,2) = norm(r_HD-rRef) / norm(rRef) * 100;
        errDat(s,i,3) = norm(r_IM-rRef) / norm(rRef) * 100;
    end
end

figure(3)
hold on
for i = 1:3
    xVar = rad2deg(taVect);
    yVar = mean(errDat(:,:,i));
    yStd =  std(errDat(:,:,i));
    plot(xVar,yVar,format{i},'MarkerFaceColor','none','LineWidth',1.25);
end
hold off
xlabel('Initial True Anomaly, $^o$')
ylabel('RMSE$(\tilde{\mathbf{r}})$, \%')
% legend('Energy','Hodograph','Improved','Location','NorthWest')
latexify(12,12,22)
setgrid
expand
% print('accuracy_ta','-dpdf','-bestfit')

% restore orbit
orbitParams = [a,e,i,o,w,f];