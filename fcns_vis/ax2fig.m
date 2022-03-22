function [xx,yy] = ax2fig(x,y)

%% Obtain arguments
hAx = gca;

%% Get limits
axun = get(hAx,'Units');
set(hAx,'Units','normalized');  % Need normaized units to do the xform
axpos = get(hAx,'Position');
axlim = axis(hAx);              % Get the axis limits [xlim ylim (zlim)]
axwidth = diff(axlim(1:2));
axheight = diff(axlim(3:4));

%% Transform data
xx = (x-axlim(1))*axpos(3)/axwidth  + axpos(1);
yy = (y-axlim(3))*axpos(4)/axheight + axpos(2);

%% Restore axes units
set(hAx,'Units',axun)

end
