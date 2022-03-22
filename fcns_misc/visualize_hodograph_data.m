function visualize_hodograph_data(vel_2d,c,R)
%VISUALIZE_HODOGRAPH_DATA Plots an orbit hodograph in 2D given its radius
%and center, then overlays the measurement data used to generate it.
%   
% Inputs:
%   vel_2d      Nx3 matrix of velocity data
%   c           hodograph center
%   R           hodograph radius

angs = linspace(0,2*pi,120)';
x = R * cos(angs) + c(1);
y = R * sin(angs) + c(2);

figure;
plot(x,y,'k')
hold on
scatter(vel_2d(:,1),vel_2d(:,2),'ro','filled')
scatter(0,0,'ko','filled')
hold off
legend('Hodograph','Data','Origin','Location','Best')
axis equal

end

