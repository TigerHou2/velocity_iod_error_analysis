%% circ2viod_obsv.m
%
% Author:
%   Tiger Hou
%
% Description:
%   Compares VIOD error trend with the error in the estimate of the
%   hodograph radius. This is the key assumption upon which the entire
%   paper is built. We should see matching trends between VIOD error and
%   hodograph radius error w.r.t. the number of observations.

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
obsvVect = 4:5:100;

mu = 1;
a = 1e5;
e = 0.9;
i = deg2rad(0);
o = deg2rad(0);
w = deg2rad(0);
f = deg2rad(0);

orbitParams = [a,e,i,o,w,f];

% total duration spanned by all measurements, as fraction of orbit period
period = 0.9;
period = period * 2*pi;
% select the nth observation's position error for comparison
selObsv = 1;
% measurement noise
noise = 3e-4;
% Monte Carlo simulation size
numSims = 3000;

% line styles
MOD = 'rx:'; % model
SIM = 'ko-.'; % simulation

% prepare measurement noise
ncube = randn(max(obsvVect),3,numSims);
ncube = ncube ./ vecnorm(ncube,2,2) .* normrnd(0,noise,1,1,numSims);

errDat = nan(numSims,length(obsvVect));
errCirc = nan(numSims,length(obsvVect));

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
    [~,R,A,B] = hodoHyp_debug(v,mu);
    
    % get position reference
    orbitParams(6) = fvect(selObsv);
    rRef = Get_Orb_Vects(orbitParams,mu);
    
ncube = randn(numObsv,3,numSims);
ncube = ncube ./ vecnorm(ncube,2,2) .* normrnd(0,noise,1,1,numSims);
    
    % Monte Carlo
    for s = 1:numSims
        nvect = ncube(1:numObsv,:,s);
        [r,Rest,Aest,Best] = hodoHyp_debug(v+nvect,mu);
        r = r(selObsv,:)';
        errDat(s,i)  = norm(r-rRef) / norm(rRef) * 100;
        
        % test svd
        [~,~,V] = svd(v+nvect,0);
        k = V(:,end);
        if k'*[0;0;1]<0
            k = -k;
        end
        
        errCirc(s,i) = 1;
        errCirc(s,i) = errCirc(s,i) * ( abs(R-Rest) / abs(R) * 100 );
%         errCirc(s,i) = errCirc(s,i) * ( norm([A-Aest,B-Best]) / abs(R) * 100 );
%         errCirc(s,i) = errCirc(s,i) * ( 1/mean(vecnorm(v+nvect,2,2)) );
%         errCirc(s,i) = errCirc(s,i) * ( 1/min(vecnorm(v+nvect,2,2)) );
%         errCirc(s,i) = errCirc(s,i) * ( norm(k-[0;0;1]) );
        
%         errDat(s,i)  = abs(R-Rest) / abs(R) * 100;
        
    end
end

xRef = obsvVect;
yRef = sqrt(mean(errDat.^2));

xVar = obsvVect;
yVar = sqrt(mean(errCirc.^2));

%% data processing & plotting
scaling = 1 / (max(yVar)-min(yVar)) * (max(yRef)-min(yRef));
yVar = yVar * scaling;
offset = - min(yVar) + min(yRef);
offset = 0;
yVar = yVar + offset;

disp(['Scaling = ' num2str(scaling)])
disp(['Offset  = ' num2str(offset)])
% disp(['ErrSpan = ' num2str(max(yRef)-min(yRef))])
% disp(['offset / errspan = ' num2str(offset/(max(yRef)-min(yRef)))])

figure;
plot(xRef,yRef,SIM,'LineWidth',1,'MarkerSize',5)
hold on
plot(xVar,yVar,MOD,'LineWidth',1,'MarkerSize',5)
hold off
y_lim = ylim .* [0.9, 1.05]; % the plot was being cut off slightly
ylim(y_lim)
% legend('RMSE$(\tilde{\mathbf{r}})$','RMSE$(\tilde{R})$','Location','Best')
xlabel('Number of Measurements')
ylabel('RMSE, \%')
latexify(10,8,16)
setgrid
expand
% svnm = [savePath 'obsvProof'];
% print(svnm,'-depsc')