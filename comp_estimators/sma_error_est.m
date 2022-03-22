%% sma_error_est.m
%
% Author:
%   Tiger Hou
%
% Description:
%   Estimates the error in VIOD position w.r.t semi-major axis length. 
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
savePath_ppt = '..\editions\space_flight_mechanics\figures_svg\';
latexify

%% setup
smaVect = linspace(1e4,1e6,15);

mu = 1;
a = 0;
e = 0.9;
i = deg2rad(0);
o = deg2rad(0);
w = deg2rad(0);
f = deg2rad(180);
    
orbitParams = [a,e,i,o,w,f];

save_plot = true;
% check if e and f0 satisfy conditions. If not, the figure is not saved.
accept(1).e = 0.1;
accept(1).ta = deg2rad(0);
accept(2).e = 0.5;
accept(2).ta = deg2rad(135);
accept(3).e = 0.9;
accept(3).ta = deg2rad(180);
if not( ( abs(accept(1).e-e) < 1e-8 && abs(accept(1).ta-f) < 1e-8 ) ...
      ||( abs(accept(2).e-e) < 1e-8 && abs(accept(2).ta-f) < 1e-8 ) ...
      ||( abs(accept(3).e-e) < 1e-8 && abs(accept(3).ta-f) < 1e-8 ) )
    save_plot = false;
    warning(['The ecc or f0 provided is not required for the ' ...
             'manuscript, therefore it will not be saved.'])
end

% total duration spanned by all measurements, as fraction of orbit period
period = 0.1;
period = period * 2*pi;
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
MOD = 'r:'; % model
SIM = 'ko-.'; % simulation

% prepare measurement noise
nGauss = normrnd(0,noise,numObsv,1,numSims);
nGauss = repmat(nGauss,1,3,1);
ncube = randn(numObsv,3,numSims);
ncube = ncube ./ vecnorm(ncube,2,2) .* nGauss;

errDat = nan(numSims,length(smaVect));
errEst = nan(1,model_res);
v = nan(numObsv,3);

%% Reference Simulation
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
    
end

% create reference sim x,y vectors
xRef = smaVect;
yRef = sqrt(mean(errDat.^2));

%% Error Model

xVar = linspace(min(smaVect),max(smaVect),model_res);

% xVar = smaVect;
% errEst = nan(1,length(xVar));

for i = 1:length(xVar)
    
    a = xVar(i);
    adj1 = sqrt(a);
    errEst(i) = adj1;
    
end

% create error model x,y vectors
yVar = errEst;

%% data processing & plotting
scaling = 1 / (max(yVar)-min(yVar)) * (max(yRef)-min(yRef));
yVar = yVar * scaling;
offset = - min(yVar) + min(yRef);
offset = 0;
yVar = yVar + offset;

disp(['Scaling = ' num2str(scaling)])
disp(['Offset  = ' num2str(offset)])

figure;
plot(xRef,yRef,SIM,'LineWidth',1.5,'MarkerSize',7)
hold on
plot(xVar,yVar,MOD,'LineWidth',1.5,'MarkerSize',7)
hold off
% legend('Monte Carlo','Error Model','Location','Best')
xlabel('Semi-Major Axis, DU')
ylabel('RMSE($\tilde{\mathbf{r}}$), \%')
latexify(10,10,16)
setgrid
expand

filename = ['smaErr_e=' num2str(e*10) '_f=' num2str(rad2deg(f))];

if save_plot
    svnm = [savePath filename];
    svnm_ppt = [savePath_ppt filename];
    print(svnm,'-depsc')
    print(svnm_ppt,'-dsvg')
end