%% rmse2df.m
%
% Author:
%   Tiger Hou
%
% Description:
%   Computes circle fitting error for different true anomaly spans.
%

%% initialize

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
e = 0;
i = deg2rad(0);
o = deg2rad(0);
w = deg2rad(0);
f = deg2rad(90);
orbitParams = [a,e,i,o,w,f];

noise = 3e-6;
numObsv = 10;
numSims = 1000;
selObsv = 1;

eVect = [0, 0.3, 0.6, 0.9];
dfVect = deg2rad(linspace(10,360,50));

ncube = randn(numObsv,3,numSims);
ncube = ncube ./ vecnorm(ncube,2,2) .* normrnd(0,noise,1,1,numSims);
v = nan(numObsv,3);

for ee = 1:length(eVect)
    
e = eVect(ee);
orbitParams(2) = e;
errDat = nan(numSims,length(dfVect));

for i = 1:length(dfVect)
    df = dfVect(i);
    f0 = f;
    f1 = f0 + df;
    E0 = 2 * atan(sqrt((1-e)/(1+e))*tan(f0/2));
    E1 = 2 * atan(sqrt((1-e)/(1+e))*tan(f1/2));
    M0 = E0 - e*sin(E0);
    M1 = E1 - e*sin(E1);
    M1 = M0 + mod(M1-M0,2*pi);
    Mvect = linspace(M0,M1,numObsv);
    Evect = kepler(Mvect,e);
    fvect = 2 * atan(sqrt((1+e)/(1-e))*tan(Evect/2));
    fvect = mod(fvect,2*pi);
    
    % ground truth data
    for j = 1:numObsv
        orbitParams(6) = fvect(j);
        [~,v(j,:)] = Get_Orb_Vects(orbitParams,mu);
    end

    % get position reference
    orbitParams(6) = fvect(selObsv);
    rRef = Get_Orb_Vects(orbitParams,mu);
    [~,R,A,B] = hodoHyp_debug(v,mu);
    
    for s = 1:numSims
        n = ncube(:,:,s);
        [~,Rest,Aest,Best] = hodoHyp_debug(v+n,mu);
        errDat(s,i) = abs(R-Rest) / abs(R) * 100;
    end
    
end

xVar = rad2deg(dfVect);
yVar = 1 ./ sqrt( sqrt(mean(errDat.^2)) );
figure(1)
hold on
plot(xVar,yVar,'x-','MarkerSize',4,'LineWidth',1.25,...
       'DisplayName',['e=' num2str(eVect(ee))])
hold off
    
end

xlabel('$\Delta{f}$, $^o$')
ylabel('RMSE($\tilde{R}$)$^{-1/2}$, \%')
legend('Location','Best')
latexify(13,10,14)
setgrid
expand
svnm = [savePath 'df'];
print(svnm,'-depsc')