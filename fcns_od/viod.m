%% viod.m
function [r] = viod(N,mu)
%VIOD Solves the three-velocity initial orbit determination problem.
%
% Author:
%   Tiger Hou
%
% Description:
%   The spacecraft's orbit is calculated from three velocity vectors
%   collected at different points in the orbit.
%
% Limitations:
%   - Assumes the measurements are provided in order.
%   - Assumes 50%+ of adjacent measurements are less than half an orbit
%   away from each other
%
% Arguments:
%   N:      [km/s]      N-by-3 matrix containing N velocity vectors
%   mu:     [km^3/s^2]  gravitational parameter of central body
%
% Sources:
%   - [1] Initial Orbit Determination from Three Velocity Vectors

%%
% preallocate matrices
%   A: [unitless]  [2*N*(N-1)]-by-[2*N] sparse matrix of linear system
%   B: [unitless]  [2*N*(N-1)]-by-1 matrix of linear system
%   C: [unitless]  [N*(N-1)/2]-by-2 matrix of combinations of
%                      two numbers selected from 1:N, non-repeating
n = size(N,1);
C = combnk(1:n,2);           % permutations of two velocity vectors i,j
comb_n = size(C,1);
% using sparse matrix to reduce memory consumption
A = spalloc(comb_n*4,n*2,7*n*(n-1)); % LHS matrix of eq. 36, see line 6
B = zeros(comb_n*4,1);               % RHS matrix of eq. 36, see line 6

% use singular value decomposition to find orbit plane normal
[~,~,V] = svd(N,0);
k = V(:,end);

% check if orbit plane normal is in the right direction
% we assume more than half of the measurements are taken less than half an
% orbit apart from its adjacent measurements.
kEst = zeros(1,3);
for i = 1:size(N,1)-1
    kEst = kEst + cross(N(i,:),N(i+1,:)) / norm(cross(N(i,:),N(i+1,:)));
end
if dot(k,kEst') < 0
    k = -k;
end

reps = size(N,1);
u = zeros(size(N')); % [3xN] matrix of u values defined in [1] eq.11
w = zeros(size(N')); % [3xN] matrix of w values defined in [1] eq.12
z = zeros(size(N')); % [3xN] matrix of z values defined in [1] eq.21
for i = 1:reps
    u(:,i) = N(i,:)' / norm(N(i,:));
    w(:,i) = cross(u(:,i),k) / norm(cross(u(:,i),k));
    z(:,i) = w(:,i) / norm(N(i,:));
end

nperms = size(C,1);
for i = 1:nperms
    % insert z and u values in top 3/4 of LHS matrix
    A((i-1)*3+1:i*3,2*C(i,1)-1) =  z(:,C(i,1));
    A((i-1)*3+1:i*3,2*C(i,1))   =  u(:,C(i,1));
    A((i-1)*3+1:i*3,2*C(i,2)-1) = -z(:,C(i,2));
    A((i-1)*3+1:i*3,2*C(i,2))   = -u(:,C(i,2));
    % insert 1 and -1 in bottom 1/4 of LHS matrix
    A((i-1)*3+nperms*3+1,2*C(i,1)-1) =  1;
    A((i-1)*3+nperms*3+1,2*C(i,2)-1) = -1;
    % insert vectors in top 3/4 of RHS matrix
    B((i-1)*3+1:i*3) = norm(N(C(i,1),:))*w(:,C(i,1)) - ...
                       norm(N(C(i,2),:))*w(:,C(i,2));
    % insert vectors in bottom 1/4 of RHS matrix
    B((i-1)*3+nperms*3+1) = norm(N(C(i,1),:))^2/2 - ...
                            norm(N(C(i,2),:))^2/2;
end
% solve the linear system
g = A\B;

a = zeros(1,reps);
b = zeros(1,reps);

a(:) = g(1:2:end);
b(:) = g(2:2:end);

r_norm = mu ./ a;

h = zeros(1,reps);
for i = 1:reps
    h(i) = sqrt( r_norm(i)^2*dot(N(i,:),N(i,:)) / ...
            (1 + dot(N(i,:),N(i,:))*(b(i)/a(i))^2) );
end
h = mean(h);

rw = zeros(size(N'));
ru = zeros(size(N'));
r  = zeros(size(N'));

for i = 1:reps
    rw(:,i) = h*z(:,i);
    ru(:,i) = b(i)*h*r_norm(i)/mu*u(:,i);
    r(:,i)  = rw(:,i) + ru(:,i);
end
r = r';
        
end %viod.m