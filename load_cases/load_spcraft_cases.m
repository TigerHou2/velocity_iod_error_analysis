%% load_spcraft_cases.m
function [noiseVect,durVect,obsVect,SMA,ECC,TA,centralBody,...
          noiseLen,durLen,obsLen] = load_spcraft_cases(varargin)
% Author:
%   Tiger Hou
%
% Description:
%   This file contains spacecraft parameter vectors that are used for case
%   studies for velocity orbit determination error. This file should be
%   called by plotSpcraftParams.m after initialization.

%% parse inputs
p = inputParser;
addOptional(p,'useGMAT',0);
parse(p,varargin{:});

%% initialize cases

% noise values
noiseVect = [1 3 5 10 15] / 1000;

% define canonical DU
DU = 1.56e5;

% measurement duration as fraction of the orbit period
durVect = [0.01 0.04 0.8]; % nd

% measurement count
obsVect = [3 10 30]; % nd

%% Earth-Neptune transfer case
mu = 1.327e11; % km^3/s^2
SMA = 2.3e9; % km
ECC = 0.5;
TA  = 0;
centralBody = 1; % Sun

%% calculate length of each parameter for GMAT
noiseLen = length(noiseVect);
durLen   = length(durVect);
obsLen   = length(obsVect);

%% calculate canonical TU
TU = sqrt(DU^3/mu);

%% convert noise to canonical units if using GMAT
if p.Results.useGMAT == 1
    noiseVect = noiseVect / DU * TU;
end

end % load_spcraft_cases.m