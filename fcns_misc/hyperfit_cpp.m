function [x,y,R] = hyperfit_cpp(dat)
%HYPERFIT_CPP Fits a circle to data using the hyperaccurate circle fit.
%   dat is an Nx2 matrix of measurements along the circle.

dat = dat';

% compute moments
means = mean(dat,2);
dat_c = dat - means;
dat_z = dat_c(1,:).^2 + dat_c(2,:).^2;
Mxx = mean(dat_c(1,:).*dat_c(1,:));
Mxy = mean(dat_c(1,:).*dat_c(2,:));
Myy = mean(dat_c(2,:).*dat_c(2,:));
Mxz = mean(dat_c(1,:).*dat_z(1,:));
Myz = mean(dat_c(2,:).*dat_z(1,:));
Mzz = mean(dat_z(1,:).*dat_z(1,:));

% compute coefficients of characteristic polynomial
Mz = Mxx + Myy;
Cov_xy = Mxx*Myy - Mxy*Mxy;
Var_z = Mzz - Mz*Mz;

A2 = 4*Cov_xy - 3*Mz*Mz - Mzz;
A1 = Var_z*Mz + 4*Cov_xy*Mz - Mxz*Mxz - Myz*Myz;
A0 = Mxz*(Mxz*Myy - Myz*Mxy) + Myz*(Myz*Mxx - Mxz*Mxy) - Var_z*Cov_xy;
A22 = A2 + A2;

% root finding using Newton's algorithm
x = 0;
y = A0;
while true
    Dy = A1 + x*(A22 + 16*x^2);
    xnew = x - y/Dy;
    if xnew == x || isinf(xnew) || isnan(xnew), break; end
    ynew = A0 + xnew*(A1 + xnew*(A2 + 4*xnew^2));
    if abs(ynew) >= abs(y), break; end
    x = xnew; y = ynew;
end

% calculate parameters of resulting circle
DET = x*x - x*Mz + Cov_xy;
Xcenter = (Mxz*(Myy - x) - Myz*Mxy)/DET/2;
Ycenter = (Myz*(Mxx - x) - Mxz*Mxy)/DET/2;

x = Xcenter + means(1);
y = Ycenter + means(2);
R = sqrt(Xcenter*Xcenter + Ycenter*Ycenter + Mz);

end

