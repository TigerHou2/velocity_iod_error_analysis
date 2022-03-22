function df = dm2df(f0,dm,e)
%DM2DF Summary of this function goes here
%   Detailed explanation goes here

E0 = 2 * atan(sqrt((1-e)/(1+e))*tan(f0/2));
M0 = E0 - e*sin(E0);
M = M0 + dm;
E = kepler(M,e);
f = 2 * atan(sqrt((1+e)/(1-e))*tan(E/2));
df = f - f0;

end

function Evect = kepler(Mvect,e)
%KEPLER solves Kepler's equation using the Conway-Laguerre algorithm.

Evect = nan(size(Mvect));

for i = 1:length(Mvect(:))
    
M = Mvect(i);

if M == 0
    Evect(i) = 0;
    continue
end

d = 1;
n = 5;

if e < 1
    E = mod(M+e/2,2*pi);
    while d > 1e-9
        E = E - n*(E-e*sin(E)-M) / ...
            (1-e*cos(E) + sign(1-e*cos(E)) ...
                        * sqrt(abs((n-1)^2*(1-e*cos(E))^2 ...
                                   -n*(n-1)*(E-e*sin(E)-M)*(e*sin(E)))));
        d = abs(E-e*sin(E)-M);
    end
else
    E = mod(M+e/2,2*pi);
    while d > 1e-9
        E = E - n*(e*sinh(E)-E-M) / ...
            (e*cosh(E)-1 + sign(e*cosh(E)-1) ...
                         * sqrt(abs((n-1)^2*(e*cosh(E)-1)^2 ...
                                    -n*(n-1)*(e*sinh(E)-E-M)*(e*sinh(E)))));
        d = abs(e*sinh(E)-E-M);
    end
end

Evect(i) = E;

end

end