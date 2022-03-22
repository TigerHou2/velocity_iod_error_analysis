%% obsv_err_est.m
%
% Author:
%   Tiger Hou
%
% Description:
%   Estimates the error in VIOD position w.r.t the number of measurements. 
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
obsvVect = 3:5:100;

mu = 1;
a = 1e5;
e = 0.5;
i = deg2rad(0);
o = deg2rad(0);
w = deg2rad(0);
f = deg2rad(90);

orbitParams = [a,e,i,o,w,f];

% total duration spanned by all measurements, as fraction of orbit period
period = 0.31;
period = period * 2*pi;
% select the nth observation's position error for comparison
selObsv = 1;
% measurement noise
noise = 3e-6;
% Monte Carlo simulation size
numSims = 3000;
% resolution for the error model
model_res = 1000;

% check if the measurement duration is 0.05, 0.3, or 0.8.
% If it is none of the above, the figure is not saved.
save_plot = true;
dur_accept = [0.05,0.3,0.8] * (2*pi);
if min(abs(period-dur_accept))>1e-8
    save_plot = false;
    warning(['The measurement duration provided is not required for ' ...
             'the manuscript, therefore it will not be saved.'])
end

% line styles
MOD = 'rx:'; % model
SIM = 'ko-.'; % simulation

% prepare measurement noise
ncube = randn(max(obsvVect),3,numSims);
ncube = ncube ./ vecnorm(ncube,2,2) * noise;

errDat = nan(numSims,length(obsvVect));
errEst = nan(1,length(obsvVect));

%% simulation
for i = 1:length(obsvVect)
    % vary number of observations
    numObsv = obsvVect(i);
    v = nan(numObsv,3);
    
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
    
    % prepare SVD error estimator
    svdErrVect = nan(numSims,1);
    
    % Monte Carlo
    for s = 1:numSims
        nvect = ncube(1:numObsv,:,s);
        r = hodoHyp(v+nvect,mu);
        r = r(selObsv,:)';
        errDat(s,i) = norm(r-rRef) / norm(rRef) * 100;
        
        % test svd
        [~,~,V] = svd(v+nvect,0);
        k = V(:,end);
        if k'*[0;0;1]<0
            k = -k;
        end
        svdErrVect(s) = norm(k-[0;0;1]);
    end
    
    % estimate error
    df = fend-f;
    df = mod(df,2*pi);
    
        % adj 1: error scales with orbit normal error
        adj1 = sqrt(mean(svdErrVect.^2));
        
        % adj 1: error scales with noise vs avg. velocity?
        adj1 = 1/(norm(vecnorm(v+nvect,2,2)))^4;
    
    errEst(i) = adj1;
    
end

%% data processing & plotting
xVar = obsvVect;
yRef = sqrt(mean(errDat.^2));
yVar = errEst;
scaling = 1 / (max(yVar)-min(yVar)) * (max(yRef)-min(yRef));
yVar = yVar * scaling;
offset = - min(yVar) + min(yRef);
yVar = yVar + offset;

disp(['Scaling = ' num2str(scaling)])
disp(['Offset  = ' num2str(offset)])

figure;
plot(xVar,yRef,SIM,'LineWidth',1,'MarkerSize',5)
hold on
plot(xVar,yVar,MOD,'LineWidth',1,'MarkerSize',5)
hold off
% legend('Simulation','Prediction','Location','Best')
xlabel('Number of Observations')
ylabel('Position MSE, \%')
latexify(10,10,16)
setgrid
expand

if save_plot
    svnm = [savePath 'obsvErr_dur=' num2str(period/2/pi*100)];
    print(svnm,'-dpdf','-bestfit')
end