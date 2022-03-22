%% hodoSuper_debug.m
function [r,R,a,b] = hodoSuper_debug(N,mu)
%HODOSUPER_DEBUG returns extra debug info compared to hodoSuper.m.
%
% Author:
%   Tiger Hou
%

% use singular value decomposition to find orbit plane normal
[~,~,V] = svd(N-mean(N),0);
k = V(:,end);

% check if orbit plane normal is in the right direction
% we assume more than half of the measurements are taken less than half an
% orbit apart from its adjacent measurements.
kEst = zeros(1,3);
for i = 1:size(N,1)-1
    kEst = kEst + cross(N(i,:),N(i+1,:));
end
if dot(k,kEst') < 0
    k = -k;
end

% define rotation matrix
ux = cross(N(1,:),k') / norm(cross(N(1,:),k')); % [1] eqn.4
uy = cross(k',ux);
T = [ux; uy; k']; % [1] eqn.5

% transform velocity from inertial frame to orbit frame
N2d = (T * N')'; % [1] eqn.3

[a,b,R] = hyperfit_cpp(N2d(:,1:2));

% find center of hodograph
c = T' * [a; b; 0]; % [1] eqn.11

% find eccentricity vector
e = cross(c,k) / R; % [1] eqn.18

% find perpendicular unit vector
upp = (N-c') ./ vecnorm(N-c',2,2); % [1] eqn.19

% find parallel unit vector
ull = cross(upp,repmat(k',size(N,1),1),2); % [1] eqn.20

% compute perpendicular velocity components [1] eqn.23
Npp = nan(size(N));
for i = 1:size(N,1)
    Npp(i,:) = ( (upp(i,:)'*upp(i,:)) * N(i,:)' )';
end

% compute position magnitudes [1] eqn.29
rho = mu * vecnorm(e'+ull,2,2) ./ vecnorm(N,2,2) ./ vecnorm(Npp,2,2);

% compute position vectors [1] eqn.30
r = rho .* ull;

end %hodoHyp_debug.m
