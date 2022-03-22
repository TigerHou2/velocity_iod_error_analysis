%% circ_span.m
%
% Author:
%   Tiger Hou
% 
% Description:
%   Models the error in the estimate for the a) Radius and b) Center of a
%   circle w.r.t. the span of measurements using the hyperfit circle 
%   fitting algorithm.
%

%% initialization

close all
clear;clc

addpath('..\fcns_misc')

R = 1;
C = [0,0];
numObsv = 17;
numSims = 1500;
noise = 5e-4;

ncube = randn(numObsv,2,numSims);
ncube = ncube ./ vecnorm(ncube,2,2) .* normrnd(0,noise,1,1,numSims);

df_range = deg2rad(linspace(5,360,200));

err_R_mat = nan(numSims,length(df_range));
err_C_mat = nan(numSims,length(df_range));

for i = 1:length(df_range)
    angs = linspace(0,df_range(i),numObsv)';
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
plot(rad2deg(df_range),1./sqrt(mean(err_R_mat)))
plot(rad2deg(df_range),1./sqrt(mean(err_C_mat)))
legend('R','C','Location','Best')