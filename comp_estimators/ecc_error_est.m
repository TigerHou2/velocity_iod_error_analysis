%% ecc_error_est.m
% 
% Author:
%   Tiger Hou
%
% Description:
%   Estimates the error in VIOD position w.r.t eccentricity. Estimated
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
eccVect = linspace(0.1,0.9,23);

mu = 1;
a = 1e5;
e = 0;
i = deg2rad(0);
o = deg2rad(0);
w = deg2rad(0);
f = deg2rad(135);

orbitParams = [a,e,i,o,w,f];

% check if initial true anomaly is 0, 135, or 180. 
% If it is none of the above, the figure is not saved.
save_plot = false;
ta_accept = deg2rad([0,135,180]);
if min(abs(f-ta_accept))>1e-8
    save_plot = false;
    warning(['The initial true anomaly provided is not required for ' ...
             'the manuscript, therefore it will not be saved.'])
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
% resolution for the error model
model_res = 1000;

% line styles
MOD = 'r:'; % model
SIM = 'ko-.'; % simulation
    
% prepare measurement noise
ncube = randn(numObsv,3,numSims);
ncube = ncube ./ vecnorm(ncube,2,2) * noise;

errDat = nan(numSims,length(eccVect));
errEst = nan(1,model_res);
v = nan(numObsv,3);

%% Reference Simulation
for i = 1:length(eccVect)
    % vary eccentricity
    e = eccVect(i);
    orbitParams(2) = e;
    
    % find measurement positions by true anomaly
    E0 = 2 * atan(sqrt((1-e)/(1+e))*tan(f/2));
    M0 = E0 - e*sin(E0);
    M = M0 + period;
    E = kepler(M,e);
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
xRef = eccVect;
yRef = sqrt(mean(errDat.^2));
y3sigma = std(errDat) * 3;
yUB = yRef + y3sigma;
yLB = yRef - y3sigma;

%% Error Model

xVar = linspace(min(eccVect),max(eccVect),model_res);
xVar = eccVect;
errEst = nan(1,length(xVar));

for i = 1:length(xVar)
    
    e = xVar(i);
    
    % adj 1: error scales inversely with hodograph radius
    adj1 = sqrt(1-e^2);
    
    % adj 2: error scales inversely with square of measurement span
    E0 = 2 * atan(sqrt((1-e)/(1+e))*tan(f/2));
    M0 = E0 - e*sin(E0);
    M = M0 + period;
    E = kepler(M,e);
    Mvect = linspace(M0,M,numObsv);
    Evect = kepler(Mvect,e);
    fvect = 2 * atan(sqrt((1+e)/(1-e))*tan(Evect/2));
    fvect = mod(fvect,2*pi);
    fend = fvect(end);
    df = fend-f;
    df = mod(df,2*pi);
    adj2 = 1 / df^2;
    
    errEst(i) = adj1 * adj2;
    
end

% create error model x,y vectors
yVar = errEst;

%% data processing & plotting
scaling = 1 / (max(yVar)-min(yVar)) * (max(yRef)-min(yRef));
yVar = yVar * scaling;
offset = - min(yVar) + min(yRef);
offset = 0;
yVar = yVar + offset;

c = 0.9 * [1,1,1];

disp(['Scaling = ' num2str(scaling)])
disp(['Offset  = ' num2str(offset)])

figure;
% patch([xRef, fliplr(xRef)], [yLB, fliplr(yUB)], c, ...
%         'EdgeColor', 'w', 'DisplayName','MC 3-$\sigma$ Bound');
hold on
% plot(xRef,yLB,'k-','DisplayName','None');
% plot(xRef,yUB,'k-','DisplayName','None');
plot(xRef,yRef,SIM,'LineWidth',1,'MarkerSize',5)
plot(xVar,yVar,MOD,'LineWidth',1.5,'MarkerSize',5)
hold off
legend('Monte Carlo','Error Model','Location','Best')
xlabel('Eccentricity')
ylabel('RMSE($\tilde{\mathbf{r}}$), \%')
latexify(10,10,16)
setgrid
expand

if save_plot
    svnm = [savePath 'eccErr_f=' num2str(rad2deg(f))];
    svnm_ppt = [savePath_ppt 'eccErr_f=' num2str(rad2deg(f))];
    print(svnm,'-depsc')
    print(svnm_ppt,'-dsvg')
end