%% dur_error_est.m
%
% Author:
%   Tiger Hou
%
% Description:
%   Estimates the error in VIOD position w.r.t measurement duration. 
%   Estimated errors are generated using analytical equations and compared
%   against Monte Carlo simulations. 

%% initialization
close all
clear;clc
addpath('..\fcns_od')
addpath('..\fcns_orb')
addpath('..\fcns_vis')
addpath('..\fcns_misc')
savePath = 'plots\';
latexify

%% setup

mu = 1;
a = 1e5;
e = 0.9;
i = deg2rad(0);
o = deg2rad(0);
w = deg2rad(0);
f = deg2rad(0);
    
orbitParams = [a,e,i,o,w,f];

% the measurement duration must be capped at df < 180 deg
% because beyond 180 deg the circle fitting error model fails.
% so we need to dynamically calculate the maximum duration allowed
% according to the eccentricity and initial true anomaly.
f0 = f;
E0 = 2 * atan(sqrt((1-e)/(1+e))*tan(f0/2));
M0 = E0 - e*sin(E0);
fmax = f0 + pi;
Emax = 2 * atan(sqrt((1-e)/(1+e))*tan(fmax/2));
Mmax = Emax - e*sin(Emax);

dM = mod(Mmax-M0,2*pi);
durVect = linspace(dM/20,dM,20);

% measurement noise
noise = 3e-6;
% number of measurements
numObsv = 10;
% Monte Carlo simulation size
numSims = 3000;
% select the nth observation's position error for comparison
selObsv = 1;
% resolution for the error model
model_res = 1000;

% line styles
MOD = 'rx:'; % model
SIM = 'ko-.'; % simulation

% prepare measurement noise
nGauss = normrnd(0,noise,numObsv,1,numSims);
nGauss = repmat(nGauss,1,3,1);
ncube = randn(numObsv,3,numSims);
ncube = ncube ./ vecnorm(ncube,2,2) .* nGauss;

errDat = nan(numSims,length(durVect));
errEst = nan(1,length(durVect));
v = nan(numObsv,3);

%% simulation
for i = 1:length(durVect)
    % vary measurement duration
    period = durVect(i);
    
    % find measurement positions by true anomaly
    E0 = 2 * atan(sqrt((1-e)/(1+e))*tan(f/2));
    M0 = E0 - e*sin(E0);
    M = M0 + period;
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
        nvect = ncube(:,:,s);
        r = hodoHyp(v+nvect,mu);
        r = r(selObsv,:)';
        errDat(s,i) = norm(r-rRef) / norm(rRef) * 100;
    end
    
    % estimate error
    df = fend-f;
    df = mod(df,2*pi);
        
        % adj 1: error scales inversely with square of measurement span
        adj1 = 1 / df^2;
    
    errEst(i) = adj1;
    
end

%% data processing & plotting
xVar = durVect /2/pi * 100;
yRef = sqrt(mean(errDat.^2));
yVar = errEst;
scaling = 1 / (max(yVar)-min(yVar)) * (max(yRef)-min(yRef));
yVar = yVar * scaling;
offset = - min(yVar) + min(yRef);
offset = 0;
yVar = yVar + offset;

disp(['Scaling = ' num2str(scaling)])
disp(['Offset  = ' num2str(offset)])

figure;
plot(xVar,yRef,SIM,'LineWidth',1,'MarkerSize',5)
hold on
plot(xVar,yVar,MOD,'LineWidth',1,'MarkerSize',5)
hold off
legend('Monte Carlo','Error Model','Location','Best')
xlabel('\% of Orbit Period')
ylabel('RMSE($\tilde{\mathbf{r}}$), \%')
latexify(10,13,18)
setgrid
expand
svnm = [savePath 'durErr_e=' num2str(e*10) '_f=' num2str(rad2deg(f))];
print(svnm,'-depsc')