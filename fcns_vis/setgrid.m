function setgrid(alpha,minoralpha,mm)
%SETGRID Sets the major and minor gridline styles for the current figure.
%   Detailed explanation goes here

if ~exist('alpha','var')
    alpha = 0.3;
end
if ~exist('minoralpha','var')
    minoralpha = 0.3;
end
if ~exist('mm','var')
    mm = 'minor';
end

grid(gca,mm)
grid on
ax = gca;
ax.GridAlpha = alpha;
ax.GridColor = [0 0 0];
ax.MinorGridAlpha = minoralpha;
ax.MinorGridColor = [0 0 0];

end

