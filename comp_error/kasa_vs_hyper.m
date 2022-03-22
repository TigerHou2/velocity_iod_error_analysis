%% kasa_vs_hyper.m
%
% Author:
%   Tiger Hou
%
% Description:
%   Compares the Kasa circle fit against the hyperaccurate circle fit and
%   generates visuals to demonstrate superiority of the hyperaccurate fit.

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

%% reference circle
R = 1;
noise = 0.01;
n_samps = 10;
c = [0,0];
plot_angs = deg2rad(1:360);
rng(1)

LW = 1.5;
SW = 1.2;

%% large arc span case
min = 0;
max = 90;
angs = deg2rad(linspace(min,max,n_samps))';
points = R * [cos(angs), sin(angs)];
n = randn(n_samps,2);
n = n ./ vecnorm(n,2,2) .* normrnd(0,noise,n_samps,1);
points = points + n;

figure(1)
plot(R*cos(plot_angs)+c(1),R*sin(plot_angs)+c(2),'r--','LineWidth',LW);
hold on
[x,y,r] = kasa(points);
plot(r*cos(plot_angs)+x,r*sin(plot_angs)+y,'b-.','LineWidth',LW);
[x,y,r] = hyperfit_cpp(points);
plot(r*cos(plot_angs)+x,r*sin(plot_angs)+y,'k','LineWidth',LW);
scatter(points(:,1),points(:,2),18,'r','LineWidth',SW);
hold off
legend('Reference','Kasa','Hyperfit','Samples',...
       'Location','Best')
xlabel('x')
ylabel('y')
axis equal
setgrid
latexify(14)

svnm = [savePath 'circfitCompare_long_captioned'];
svnm_ppt = [savePath_ppt 'circfitCompare_long'];
print(svnm,'-depsc')
print(svnm_ppt,'-dsvg')

%% small arc span case
min = 0;
max = 25;
angs = deg2rad(linspace(min,max,n_samps))';
points = R * [cos(angs), sin(angs)];
n = randn(n_samps,2);
n = n ./ vecnorm(n,2,2) .* normrnd(0,noise,n_samps,1);
points = points + n;

figure(2)
plot(R*cos(plot_angs)+c(1),R*sin(plot_angs)+c(2),'r--','LineWidth',LW);
hold on
[x,y,r] = kasa(points);
plot(r*cos(plot_angs)+x,r*sin(plot_angs)+y,'b-.','LineWidth',LW);
[x,y,r] = hyperfit_cpp(points);
plot(r*cos(plot_angs)+x,r*sin(plot_angs)+y,'k','LineWidth',LW);
scatter(points(:,1),points(:,2),18,'r','LineWidth',SW);
hold off
% legend('Reference','Kasa','Hyperfit','Samples',...
%        'Location','NorthEastOutside')
xlabel('x')
ylabel('y')
axis equal
setgrid
latexify(14)

svnm = [savePath 'circfitCompare_short'];
svnm_ppt = [savePath_ppt 'circfitCompare_short'];
print(svnm,'-depsc')
print(svnm_ppt,'-dsvg')