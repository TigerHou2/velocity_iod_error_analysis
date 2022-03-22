%% hodoAvg.m
function [r] = hodoAvg(vel,time,mu)
%HODOHYP Solves the three-velocity initial orbit determination problem.
%
% Author:
%   Tiger Hou
%
% Note:
%   This version uses John Christian's StarNAV paper:
%   https://doi.org/10.3390/s19194064
%   This version also applies the hyperaccurate circle fit from the
%   previous version of the algorithm, hodoHyp.m
%
% Description:
%   The spacecraft's orbit is calculated using the velocity hodograph from
%   three velocity vectors collected at different points in orbit.
%
% Limitations:
%   - Assumes the measurements are provided in order.
%   - Assumes 50%+ of adjacent measurements are less than half an orbit
%   away from each other
%
% Arguments:
%   vel:    [km/s]      N-by-3 matrix containing N velocity vectors
%   time:   [s]         N-by-1 matrix containing relative measurement times
%   mu:     [km^3/s^2]  gravitational parameter of central body
%
% References:
%   [1] - Geometric Solutions for Problems in Velocity-Based Orbit Determination
%   [2] - StarNAV paper: https://doi.org/10.3390/s19194064

% transpose N to make vectors columns
Vel = vel';

% use singular value decomposition to find orbit plane
[~,~,V] = svd(Vel',0);
k = V(:,end);

% count the number of velocity measurements
n = size(Vel,2);

% check if orbit plane normal is in the right direction
% we assume more than half of the measurements are taken less than half an
% orbit apart from their adjacent measurements.
kEst = zeros(3,1);
for i = 1:n-1
    kEst = kEst + cross(Vel(:,i),Vel(:,i+1));
end
if (k'*kEst) < 0
    k = -k;
end

% define rotation matrix from inertial frame to orbit frame
ux = cross(Vel(:,1),k) / norm(cross(Vel(:,1),k)); % [1] eqn.4
uy = cross(k,ux);
T = [ux; uy; k]; % [1] eqn.5

% transform velocity from inertial frame to orbit frame
Vel_orb = T * Vel; % [1] eqn.3

% perform circle fit to find the hodograph in the orbit frame
[a,b,R] = hyperfit_cpp(Vel_orb(:,1:2));

% find center of hodograph in the inertial frame
c = T' * [a; b; 0]; % [1] eqn.11

% find eccentricity vector in the inertial frame
e = cross(c,k) / R; % [1] eqn.18

% find the semi-major axis
% StarNAV eqn. 149
a = mu / (R^2-c'*c);

% normalize the hodograp center and velcity vectors
Vnorm = Vel ./ vecnorm(Vel,2,2);
cnorm = c / norm(c);

% create the cross product of c and N (velocity)
C = repmat(c,1,n); % Nx3 matrix
Cn = repmat(cnorm,1,n); % Nx3 matrix
cross_Cn_Nn = cross(Cn,Vnorm,1); % Nx3 matrix
cross_C_Nn = cross(C,Vnorm,1); % Nx3 matrix

% find the angle phi between the hodograph center and the velocity
% StarNAV eqn. 134
phi = acos( cnorm' * Vnorm); % 1xN matrix

% find the corresponding true anomalies
% StarNAV eqn. 137
fval = phi+asin(vecnorm(cross_C_Nn/R,2,1)); % 1xN matrix
f = double( k'*cross_Cn_Nn> 0 ).* fval ...
  + double( k'*cross_Cn_Nn<=0 ).* (2*pi - fval); % 1xN matrix

% find the mean time of periapsis passage
sin_E = sin(f)*sqrt(1-norm(e)^2) / (1+norm(e)*cos(f)); % eqn. 138
cos_E = (norm(e)+cos(f)) / (1+norm(e)*cos(f));         % eqn. 139
E = atan2(sin_E,cos_E);
M = E - norm(e)*sin_E;                                 % eqn. 140
tpo = mean(time'-sqrt(a^3/mu)*M); % 1xN matrix            % eqn. 142

% find the averaged mean anomaly
% StarNAV eqn. 143
Mavg = sqrt(mu/a^3)*(time'-tpo);
Eavg = kepler(Mavg,norm(e));

% find the averaged true anomaly
sin_f = sin(Eavg)*sqrt(1-norm(e)^2) / (1-norm(e)*cos(Eavg)); % eqn. 144
cos_f = (cos(Eavg)-norm(e)) / (1-norm(e)*cos(Eavg));         % eqn. 145

% find the unit position vectors
% StarNAV eqn. 146
rnorm = cos_f*e + sin_f*cross(k,e/norm(e));

% find the position magnitudes
% StarNAV eqn. 147
rho = mu/R^2 / (1+norm(e)*cos_f);

% find the position vectors
r = rnorm .* rho;

end %hodoAvg.m
