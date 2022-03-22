%% hodoHyp_debug.m
function [r,R,a,b,N2d,vel2d] = hodoHyp_debug(v_noisy,mu,v_truth)
%HODOHYP_DEBUG returns extra debug info compared to hodoHyp.m.
%
% Author:
%   Tiger Hou
%

%=====================================================
%DEBUG
if (nargin==2)
    N1 = v_noisy;
    N2 = v_noisy;
else
    N1 = v_noisy;  % circle fitting phase
    N2 = v_noisy;  % position calc phase
end

%=====================================================

% use singular value decomposition to find orbit plane normal
[~,~,V] = svd(N1,0);
k = V(:,end);

% check if orbit plane normal is in the right direction
% we assume more than half of the measurements are taken less than half an
% orbit apart from its adjacent measurements.
kEst = zeros(1,3);
for i = 1:size(N1,1)-1
    kEst = kEst + cross(N1(i,:),N1(i+1,:));
end
if dot(k,kEst') < 0
    k = -k;
end

% define rotation matrix
ux = cross(N1(1,:),k') / norm(cross(N1(1,:),k')); % [1] eqn.4
uy = cross(k',ux);
T = [ux; uy; k']; % [1] eqn.5

% transform velocity from inertial frame to orbit frame
N2d = (T * N1')'; % [1] eqn.3

[a,b,R] = hyperfit(N2d(:,1:2));

% find center of hodograph
c_2d = [a,b,0]';
c = T' * c_2d; % [1] eqn.11

% find eccentricity vector
e = cross(c,k) / R; % [1] eqn.18

%=====================================================
%TESTING
% instead of using the actual velocity data we measured to find the
% position vectors, we will use the point nearest on the predicted
% hodograph instead.

if 1

N2d = (T * N2')';
direction = N2d - c_2d';
vel2d = R * direction ./ vecnorm(direction,2,2) + c_2d';
N2 = (T' * vel2d')';

r = nan(size(N2));
% BELOW are three forms for solving r that yield similar results

% for i = 1:size(N2,1)
%     r(i,:) = mu/R ...
%            * sin(pi/2-atan2(norm(cross(N2(i,:),N1(i,:))),...
%                                    dot(N2(i,:),N1(i,:)) ) ) ...
%            / ( N2(i,:) * (N2(i,:)-c')' ) ...
%            * cross(N2(i,:)-c',k);
% end
% for i = 1:size(N2,1)
%     f = atan2(norm(cross(N2d(i,:)-[a,b,0],[a,b,0])),...
%                      dot(N2d(i,:)-[a,b,0],[a,b,0]));
%     upp = (N2(i,:)-c') / norm(N2(i,:)-c');
%     r(i,:) = mu/R ...
%            / sqrt(norm(N2(i,:))^2-norm(c)^2*sin(f)^2) ...
%            * cross(upp,k');
% end
for i = 1:size(N2,1)
    vpp = (N2(i,:)-c');
    r(i,:) = mu/R * cross(vpp,k') / (N2(i,:)*vpp');
end

return

end


%=====================================================

% find perpendicular unit vector
upp = (N2-c') ./ vecnorm(N2-c',2,2); % [1] eqn.19

% find parallel unit vector
ull = cross(upp,repmat(k',size(N2,1),1),2); % [1] eqn.20

% compute perpendicular velocity components [1] eqn.23
Npp = nan(size(N2));
for i = 1:size(N2,1)
    Npp(i,:) = ( (upp(i,:)'*upp(i,:)) * N2(i,:)' )';
end

% compute position magnitudes [1] eqn.29
rho = mu * vecnorm(e'+ull,2,2) ./ vecnorm(N2,2,2) ./ vecnorm(Npp,2,2);

% compute position vectors [1] eqn.30
r = rho .* ull;

end %hodoHyp_debug.m
