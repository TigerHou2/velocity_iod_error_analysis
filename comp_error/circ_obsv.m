%% circ_obsv.m
%
% Author:
%   Tiger Hou
% 
% Description:
%   Models the error in the estimate for the a) Radius and b) Center of a
%   circle w.r.t. the number of observations using the hyperfit circle 
%   fitting algorithm.
%

%% initialization

close all
clear;clc

addpath('..\fcns_misc')

R = 1;
C = [0,0];
df = deg2rad(90);
numSims = 1500;
noise = 5e-4;

obsvVect = 3:5:80;

err_R_mat = nan(numSims,length(obsvVect));
err_C_mat = nan(numSims,length(obsvVect));

for i = 1:length(obsvVect)
    
    numObsv = obsvVect(i);
    
    ncube = randn(numObsv,2,numSims);
    ncube = ncube ./ vecnorm(ncube,2,2) .* normrnd(0,noise,1,1,numSims);
    
    angs = linspace(0,df,numObsv)';
    xy = [cos(angs), sin(angs)] * R + C;
    
    for s = 1:numSims
        n = ncube(:,:,s);
        [a,b,Rest] = hyperfit_cpp(xy+n);
        err_R_mat(s,i) = abs(R-Rest) / abs(R) * 100;
        err_C_mat(s,i) = norm([a,b]-C) / abs(R) * 100;
    end
    
end

figure(1)
hold on
plot(obsvVect, 1./sqrt(mean(err_R_mat)))
plot(obsvVect, 1./sqrt(mean(err_C_mat)))
legend('R','C','Location','Best')