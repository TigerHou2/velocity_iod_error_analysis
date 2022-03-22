%% hodoHyp.m
function [r,R] = hodo_timed(vel,t,mu)
%HODOHYP Solves the three-velocity initial orbit determination problem.
%
% Author:
%   Tiger Hou
%
% Note:
%   This version continues from hodoHyp.m and adds timing data to each
%   measurement, using John Christian's StarNAV paper to find the mean time
%   of periapsis passage. 
%
% Arguments:
%   N:      [km/s]      N-by-3 matrix containing N velocity vectors
%   t:      [s]         N-by-1 vector of measurement times for velocities
%   mu:     [km^3/s^2]  gravitational parameter of central body
%
% References:
%   [1] - Geometric Solutions for Problems in Velocity-Based Orbit Determination
%   [2] - StarNAV paper

% use singular value decomposition to find orbit plane normal
[~,~,V] = svd(vel,0);
k = V(:,end);

% check if orbit plane normal is in the right direction
% we assume more than half of the measurements are taken less than half an
% orbit apart from its adjacent measurements.
kEst = zeros(1,3);
for i = 1:size(vel,1)-1
    kEst = kEst + cross(vel(i,:),vel(i+1,:));
end
if dot(k,kEst') < 0
    k = -k;
end

% define rotation matrix
ux = cross(vel(1,:),k') / norm(cross(vel(1,:),k')); % [1] eqn.4
uy = cross(k',ux);
T = [ux; uy; k']; % [1] eqn.5

% transform velocity from inertial frame to orbit frame
vel2d = (T * vel')'; % [1] eqn.3

[a,b,R] = hyperfit_cpp(vel2d(:,1:2));

% find center of hodograph
c = T' * [a; b; 0]; % [1] eqn.11

% find eccentricity vector
e = cross(c,k) / R; % [1] eqn.18

% find phi, the angle between the hodograph center and current velocity
uc = c / norm(c);
uv = vel ./ vecnorm(vel,2,2);
phi = acos(uv*uc); % [2] eqn.136

% temporary patch for e < 0 or e >= 1
if norm(e) < 0
    e = e / norm(e) * 1e-3;
elseif norm(e) >= 1
    e = e / norm(e) * 0.999;
end

% find true anomalies
f = nan(size(phi));
for i = 1:length(f)
    if cross(uc,uv(i,:))*k > 0  % [2] eqn.137
        f(i) = phi(i) + asin(1/R*norm(cross(c,uv(i,:))));
    else
        f(i) = 2*pi - phi(i) - asin(1/R*norm(cross(c,uv(i,:))));
    end
end

% find mean anomaly
E = 2*atan( tan(f/2) * sqrt((1-norm(e))/(1+norm(e))) );
M = E - norm(e)*sin(E); % [2] eqn.138-140
M = mod(M,2*pi);

% find semi-major axis
a = mu / R^2 / (1-norm(e)^2);

% find mean time of periapsis passage
tp0 = 1/size(vel,1) * sum(t-sqrt(a^3/mu)*M); % [2] eqn.142

% find mean anomaly since periapsis passage for all measurements
Mm = sqrt(mu/a^3) * (t-tp0); % [2] eqn.143

% find true anoamly since periapsis passage for all measurements
Em = kepler(Mm,norm(e));
Fm = 2*atan( tan(Em/2) * sqrt((1+norm(e))/(1-norm(e))) );

% find position unit vectors
ur = cos(Fm)*e'/norm(e) + sin(Fm)*cross(k,e/norm(e))'; % [2] eqn.146

% find position magnitudes
rho = mu/R^2./(1+norm(e)*cos(Fm)); % [2] eqn.147

% find position vectors
r = rho .* ur; % [2] eqn.130

end %hodo_timed.m
