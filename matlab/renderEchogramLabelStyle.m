function h = renderEchogramLabelStyle(ax, x, depth, NASC)
%RENDERECHOGRAMLABELSTYLE Draw an echogram the way the labeling tool shows it:
% linear NASC, color range capped at the 98th percentile of positive values
% (so a few extreme cells do not wash out the structure), and the cmocean
% 'deep' colormap if cmocean is on the path (parula otherwise).
%
% INPUTS
%   ax     target axes
%   x      [1 x nCol] or [nCol x 1] time/index coordinate
%   depth  [1 x nDepth] or [nDepth x 1] depth coordinate, ascending
%   NASC   [nCol x nDepth] linear NASC
%
% OUTPUT  h, the pcolor surface handle.

h = pcolorLikeImagesc(ax, x, depth, NASC');
set(ax, 'YDir', 'reverse');
if exist('cmocean', 'file') == 2
    colormap(ax, flipud(cmocean('deep', 256)));
else
    colormap(ax, parula(256));
end
finPos = NASC(isfinite(NASC) & NASC > 0);
if ~isempty(finPos)
    clim(ax, [0, prctile(finPos, 98)]);
end
end
