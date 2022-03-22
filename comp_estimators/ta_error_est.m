%% ta_error_est.m
%
% Author:
%   Tiger Hou
%
% Description:
%   Estimates the error in VIOD position w.r.t true anomaly. Estimated
%   errors are generated using analytical equations and compared against
%   Monte Carlo simulations. 

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
taVect = deg2rad(linspace(0,360,37));
taVect(end) = [];

mu = 1;
a = 1e5;
e = 0.9;
i = deg2rad(0);
o = deg2rad(0);
w = deg2rad(0);
f = deg2rad(0);

orbitParams = [a,e,i,o,w,f];

% check if eccentricity is 0.1, 0.5, or 0.9. If it is none of the above,
% the figure is not saved.
save_plot = false;
ecc_accept = [0.1,0.5,0.9];
if min(abs(e-ecc_accept))>1e-8
    save_plot = false;
    warning(['The eccentricity provided is not required for the ' ...
             'manuscript, therefore it will not be saved.'])
end

% total duration spanned by all measurements, as fraction of orbit period
period = 0.1;
period = period * 2*pi;
% select the nth observation's position error for comparison
selObsv = 1;
% number of measurements
numObsv = 10;
% measurement noise
noise = 3e-6;
% Monte Carlo simulation size
numSims = 3000;

% line styles
MOD = 'r:'; % model
SIM = 'ko-.'; % simulation
    
% prepare measurement noise
ncube = randn(numObsv,3,numSims);
ncube = ncube ./ vecnorm(ncube,2,2) * noise;

errDat = nan(numSims,length(taVect));
errEst = nan(1,length(taVect));
v = nan(numObsv,3);
 
%% Reference Simulation
for i = 1:length(taVect)
    % vary true anomaly
    f0 = taVect(i);
    orbitParams(6) = f0;
    
    % find measurement positions by true anomaly
    E0 = 2 * atan(sqrt((1-e)/(1+e))*tan(f0/2));
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

% create reference simulation x,y vectors
xRef = taVect;
yRef = sqrt(mean(errDat.^2));

%% Error Model

model_resolution = 1000;
xVar = linspace(min(taVect),max(taVect),model_resolution);

xVar = taVect;
errEst = nan(1,length(xVar));

for i = 1:length(xVar)
    
    f0 = xVar(i);
    
    % adj 1: error scales inversely with square of measurement span
    E0 = 2 * atan(sqrt((1-e)/(1+e))*tan(f0/2));
    M0 = E0 - e*sin(E0);
    Me = M0 + period;
    Ee = kepler(Me,e);
    fe = 2 * atan(sqrt((1+e)/(1-e))*tan(Ee/2));
    fe = mod(fe,2*pi);
    df = fe-f0;
    df = mod(df,2*pi);
    adj1 = 1 / df^2;
    
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
plot(rad2deg(xRef),yRef,SIM,'LineWidth',1,'MarkerSize',5)
hold on
plot(rad2deg(xVar),yVar,MOD,'LineWidth',1.5,'MarkerSize',5)
hold off
xlim([0,360])
% legend('Monte Carlo','Error Model','Location','Best')
xlabel('Initial True Anomaly, $^o$')
ylabel('RMSE($\tilde{\mathbf{r}}$), \%')
latexify(10,10,16)
setgrid
expand

if save_plot
    svnm = [savePath 'taErr_e=' num2str(e*10)];
    svnm_ppt = [savePath_ppt 'taErr_e=' num2str(e*10)];
    print(svnm,'-depsc')
    print(svnm_ppt,'-dsvg')
end