function [a,e,i,omg,w,f] = Get_Orb_Params(r,v,mu)
%GET_ORB_PARAMS Summary of this function goes here
%   Detailed explanation goes here

a = norm(r)/(2 - norm(r)*dot(v,v)/mu);

e = ((dot(v,v)-mu/norm(r))*r - dot(r,v)*v)/mu;

I = [1 0 0];
J = [0 1 0];
K = [0 0 1];

h = cross(r,v);

n = cross(K, h/norm(h));

i = acos(dot(h,K)/norm(h));

omg = acos(dot(n,I)/norm(n));
if dot(n,J) < 0
    omg = 2*pi - omg;
end

w = acos(dot(n,e)/norm(n)/norm(e));
if dot(e,K) < 0
    w = 2*pi - w;
end
   
f = acos(dot(e,r)/norm(e)/norm(r));
if dot(r,v) < 0
    f = 2*pi - f;
end

end

