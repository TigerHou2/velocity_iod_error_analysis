%% load_cases_orbit.m
function [smaVect,eccVect,taVect,noise_canon,names,centralBody,...
          smaLen, eccLen, taLen, noise, numObsv] = load_orbit_cases(varargin)
% Author:
%   Tiger Hou
%
% Description:
%   This file contains orbital parameter vectors that are used for case
%   studies for velocity orbit determination error. This file should be
%   called by plotOrbitParams.m after initialization.

%% parse inputs
p = inputParser;
addOptional(p,'useGMAT',0);
addOptional(p,'eccVect',[0.1,0.5,0.9]);
addOptional(p,'taVect' ,[0, 20, 45, 90, 120, 160, 180]);
parse(p,varargin{:});

%% define noise, initialize cases

% we wish to set up canonical units such that
% noise is 1e-6 DU/TU
noise_canon = 1e-6;

noise = 3; % m/s
noise = noise / 1000; % km/s

eccVect = p.Results.eccVect;
taVect = p.Results.taVect;
taVect  = deg2rad(taVect);
smaVect = [];
names = {};
centralBody = [];

%% GEO case
mu  = 3.986e5; % km^3/s^2
SMA = 41864; % km
if p.Results.useGMAT == 1
    smaVect(end+1) = SMA;
    names = 0;
else
    smaVect(end+1) = SMA / (mu/noise^2) / noise_canon^2;
    names{end+1} = 'GEO';
end
centralBody(end+1) = 0;

%% Earth-Neptune transfer case
mu = 1.327e11; % km^3/s^2
SMA = 2.3e9; % km
if p.Results.useGMAT == 1
    smaVect(end+1) = SMA;
    names = 0;
else
    smaVect(end+1) = SMA / (mu/noise^2) / noise_canon^2;
    names{end+1} = 'Earth-Neptune';
end
centralBody(end+1) = 1;

%% Earth-Mars transfer case
mu = 1.327e11; % km^3/s^2
SMA = 1.9e8; % km
if p.Results.useGMAT == 1
    smaVect(end+1) = SMA;
    names = 0;
else
    smaVect(end+1) = SMA / (mu/noise^2) / noise_canon^2;
    names{end+1} = 'Earth-Mars';
end
centralBody(end+1) = 1;

%% calculate length of each parameter for GMAT
smaLen = length(smaVect);
eccLen = length(eccVect);
taLen  = length(taVect);

%% number of observations
numObsv = 3;

end % load_orbit_cases.m