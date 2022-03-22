%% circ2viod_dur.m
%
% Author:
%   Tiger Hou
%
% Description:
%   Compares VIOD error trend with the error in the estimate of the
%   hodograph radius. This is the key assumption upon which the entire
%   paper is built. We should see matching trends between VIOD error and
%   hodograph radius error w.r.t. the total measurement time interval.

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
durVect = linspace(0.1,0.9,25);
durVect = durVect * 2*pi;

mu = 1;
a = 1e5;
e = 0.8;
i = deg2rad(0);
o = deg2rad(0);
w = deg2rad(0);
f = deg2rad(0);
    
orbitParams = [a,e,i,o,w,f];

% measurement noise
noise = 8e-4;
% number of measurements
numObsv = 10;
% Monte Carlo simulation size
numSims = 3000;
% select the nth observation's position error for comparison
selObsv = 1;

% line styles
MOD = 'rx:'; % model
SIM = 'ko-.'; % simulation

% prepare measurement noise
nGauss = normrnd(0,noise,numObsv,1,numSims);
nGauss = repmat(nGauss,1,3,1);
ncube = randn(numObsv,3,numSims);
ncube = ncube ./ vecnorm(ncube,2,2) .* nGauss;

errDat = nan(numSims,length(durVect));
errCirc = nan(numSims,length(durVect));
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
    [~,R,a,b,vel] = hodoHyp_debug(v,mu);
    
    % get position reference
    orbitParams(6) = fvect(selObsv);
    rRef = Get_Orb_Vects(orbitParams,mu);
    
    % Monte Carlo
    for s = 1:numSims
        nvect = ncube(:,:,s);
        [r,Rest,a_est,b_est,vel_2d] = hodoHyp_debug(v+nvect,mu);
        r = r(selObsv,:)';
        
        errDat(s,i)  = norm(r-rRef) / norm(rRef) * 100;
        errCirc(s,i) = abs(R-Rest) / abs(R) * 100;
%         errCirc(s,i) = norm([a_est,b_est]-[a,b]) / norm([a,b]) * 100;
        
    end
end

xRef = durVect / (2*pi) * 100;
yRef = sqrt(mean(errDat.^2));

xVar = durVect / (2*pi) * 100;
yVar = sqrt(mean(errCirc.^2));

%% data processing & plotting
scaling = 1 / (max(yVar)-min(yVar)) * (max(yRef)-min(yRef));
scaling = 1;
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
y_lim = ylim .* [0.9, 1.1]; % the plot was being cut off slightly
ylim(y_lim)
legend('RMSE$(\tilde{\mathbf{r}})$','RMSE$(\tilde{R})$','Location','Best')
xlabel('\% of Orbit Period')
ylabel('RMSE, \%')
latexify(10,8,16)
setgrid
expand
% svnm = [savePath 'durProof'];
% print(svnm,'-depsc')