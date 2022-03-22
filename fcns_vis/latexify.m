function latexify(v1,v2,v3)
%LATEXIFY Sets text interpreter to latex for all figures.
%   Optional arguments also allow custom font and figure sizing.
%   Options:
%       # args = 1: sets font size (v1)
%       # args = 2: sets figure width (v1) and height (v2)
%       # args = 3: sets figure width (v1), height (v2), and font size (v3)
%
% Author:
%   Tiger Hou

%%
set(groot,'defaulttextinterpreter','latex');  
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

g = groot;
if ~isempty(g.Children) % only set figure properties if figure exists
    set(gca, 'Color', [1 1 1])
    set(gcf, 'PaperPositionMode', 'manual')
    set(gcf, 'Color', [1 1 1])
    set(gcf, 'PaperUnits', 'centimeters')
    set(gcf, 'Units', 'centimeters' )

    switch nargin
        case 1
            fontSize = v1;
            set(findall(gcf,'-property','FontSize'),'FontSize',fontSize)
        case 2
            w = v1;
            h = v2;
            set(gcf, 'PaperSize', [w h])
            set(gcf, 'Position', [0.2 1.2 w h])
            set(gcf, 'PaperPosition', [0.2 1.2 w h])
        case 3
            w = v1;
            h = v2;
            fontSize = v3;
            set(gcf, 'PaperSize', [w h])
            set(gcf, 'Position', [0.2 1.2 w h])
            set(gcf, 'PaperPosition', [0.2 1.2 w h])
            set(findall(gcf,'-property','FontSize'),'FontSize',fontSize)
    end
end

end %latexify.m
