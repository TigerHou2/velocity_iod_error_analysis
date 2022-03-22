function [x,y,R] = kasa(dat)
%KASA Fits a circle to data using the Kasa circle fit.
%   dat is an Nx2 matrix of measurements along the circle.

% construct linear system
A = 2*dat;
A(:,3) = -1;
B = dat.^2;
B = sum(B(:,1:2),2);
x = A\B;

% find radius of hodograph
R = sqrt(x(1)^2 + x(2)^2 - x(3));

y = x(2);
x = x(1);

end

