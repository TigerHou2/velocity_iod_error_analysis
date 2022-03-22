%% circ2sma.m
%
% Author:
%   Tiger Hou
%
% Description:
%   Compares VIOD error trend with the error in the estimate of the
%   hodograph radius. This is the key assumption upon which the entire
%   paper is built. We should see matching trends between VIOD error and
%   hodograph radius error w.r.t. semi-major axis.

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

rng(5)

smaVect = linspace(1e4,1e6,15);

mu = 1;
a = 0;
e = 0.9;
i = deg2rad(0);
o = deg2rad(0);
w = deg2rad(0);
f = deg2rad(190);
    
orbitParams = [a,e,i,o,w,f];

% total duration spanned by all measurements, as fraction of orbit period
period = 0.9;
period = period * 2*pi;
% measurement noise
noise = 12e-5;
% number of measurements
numObsv = 10;
% Monte Carlo simulation size
numSims = 1000;
% select the nth observation's position error for comparison
selObsv = 1;

% line styles
MOD = 'rx:'; % model
SIM = 'ko-.'; % simulation

% prepare measurement noise
ncube = randn(numObsv,3,numSims);
ncube = ncube ./ vecnorm(ncube,2,2) .* normrnd(0,noise,1,1,numSims);

errDat = nan(numSims,length(smaVect));
errCirc = nan(numSims,length(smaVect));
v = nan(numObsv,3);

%% simulation
for i = 1:length(smaVect)
    % vary semi-major axis length
    sma = smaVect(i);
    orbitParams(1) = sma;
    
    % find measurement positions by true anomaly
    E0 = 2 * atan(sqrt((1-e)/(1+e))*tan(f/2));
    M0 = E0 - e*sin(E0);
    M = M0 + period;
    Mvect = linspace(M0,M,numObsv);
    Evect = kepler(Mvect,e);
    fvect = 2 * atan(sqrt((1+e)/(1-e))*tan(Evect/2));
    fvect = mod(fvect,2*pi);
    fend = fvect(end);
    
    if i==1
        df = rad2deg(mod(fend-f,2*pi));
        disp(['df = ' num2str(df) ' deg']);
    end

    % ground truth data
    for j = 1:numObsv
        orbitParams(6) = fvect(j);
        [~,v(j,:)] = Get_Orb_Vects(orbitParams,mu);
    end
    
    [~,R,a,b,vel] = hodoHyp_debug(v,mu);

    % get position reference
    orbitParams(6) = fvect(selObsv);
    rRef = Get_Orb_Vects(orbitParams,mu);
    
    % Monte Carlo
    for s = 1:numSims
        nvect = ncube(:,:,s);
        [r,Rest,a_est,b_est,dat2d,datProj] = hodoHyp_debug(v+nvect,mu);
        r = r(selObsv,:)';
        
        if (Rest > R*4) || (Rest < R/4)
            continue
        end
        
        errDat(s,i) = norm(r-rRef) / norm(rRef) * 100;
        errCirc(s,i) = abs(R-Rest) / abs(R) * 100;
        
    end
    
end

xRef = smaVect;
yRef = sqrt(mean(errDat.^2,'omitnan'));

xVar = smaVect;
yVar = sqrt(mean(errCirc.^2,'omitnan'));

if max(sum(isnan(errDat))) > 0
    disp([  'Removed ' num2str(max(sum(isnan(errDat)))) ...
            ' extreme cases with large errors.'])
end

%% data processing & plotting
scaling = 1 / (max(yVar)-min(yVar)) * (max(yRef)-min(yRef));
yVar = yVar * scaling;
offset = - min(yVar) + min(yRef);
offset = 0;
yVar = yVar + offset;

disp(['Scaling = ' num2str(scaling)])
disp(['Offset  = ' num2str(offset)])

figure;
plot(xRef,yRef,SIM,'LineWidth',1,'MarkerSize',5)
hold on
plot(xVar,yVar,MOD,'LineWidth',1,'MarkerSize',5)
hold off
legend('RMSE$(\tilde{\mathbf{r}})$','RMSE$(\tilde{R})$','Location','Best')
xlabel('Semi-Major Axis, DU')
ylabel('RMSE, \%')
latexify(10,8,15)
setgrid
expand
% svnm = [savePath 'smaProof'];
% print(svnm,'-depsc')