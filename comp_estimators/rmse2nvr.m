%% rmse2nvr.m
%
% Author:
%   Tiger Hou
%
% Description:
%   Computes circle fitting error for different velocity-to-noise ratios.
%

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
e = 0;
i = deg2rad(0);
o = deg2rad(0);
w = deg2rad(0);
f = deg2rad(0);
orbitParams = [a,e,i,o,w,f];

dM = 0.1 * (2*pi);
numObsv = 10;
numSims = 1000;
selObsv = 1;

nvrVect = linspace(1e-7,1e-5,30);
ncube = randn(numObsv,3,numSims);
ncube = ncube ./ vecnorm(ncube,2,2) .* normrnd(0,1,1,1,numSims);
v = nan(numObsv,3);

eVect = [0.3 0.8];
fVect = deg2rad([0 90 135]);

for ee = 1:length(eVect)
for ff = 1:length(fVect)
    
if ee == 2 && ff == 3
    continue
end
    
e = eVect(ee);
f0 = fVect(ff);
orbitParams(2) = e;
orbitParams(6) = f0;

% find measurement positions by true anomaly
E0 = 2 * atan(sqrt((1-e)/(1+e))*tan(f0/2));
M0 = E0 - e*sin(E0);
M = M0 + dM;
Mvect = linspace(M0,M,numObsv);
Evect = kepler(Mvect,e);
fvect = 2 * atan(sqrt((1+e)/(1-e))*tan(Evect/2));
fvect = mod(fvect,2*pi);
fend = fvect(end);

% ground truth data
for i = 1:numObsv
    orbitParams(6) = fvect(i);
    [~,v(i,:)] = Get_Orb_Vects(orbitParams,mu);
end

% get position reference
orbitParams(6) = fvect(selObsv);
rRef = Get_Orb_Vects(orbitParams,mu);
[~,R,A,B] = hodoHyp_debug(v,mu);

errDat = nan(numSims,length(nvrVect));

for i = 1:length(nvrVect)
    noise = nvrVect(i);
    ncube_loc = ncube * noise;
    for s = 1:numSims
        n = ncube_loc(:,:,s);
        [~,Rest,Aest,Best] = hodoHyp_debug(v+n,mu);
        errDat(s,i) = abs(R-Rest) / abs(R) * 100;
    end
end

xVar = nvrVect;
yVar = sqrt(mean(errDat.^2));
figure(1)
hold on
plot(xVar,yVar,'x-','MarkerSize',4,'LineWidth',1.25,...
       'DisplayName',['e=' num2str(e) ', f0=' num2str(rad2deg(f0)) '$^o$'])
hold off

end
end

xlabel('$\sigma/R$, $TU^{-1}$')
ylabel('RMSE($\tilde{R}$), \%')
legend('Location','Best')
latexify(13,10,14)
setgrid
expand
svnm = [savePath 'NVR'];
print(svnm,'-depsc')