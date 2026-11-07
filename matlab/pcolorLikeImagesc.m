function h = pcolorLikeImagesc(ax, x, y, C)
%PCOLORLIKEIMAGESC Draw C like imagesc(ax,x,y,C) would (x along columns, y along
% rows, every cell shown), but as a pcolor surface so that vector export
% (SVG/EPS) keeps one path per cell instead of a single raster image.
%
% pcolor treats x and y as cell-boundary vertices and would drop the last row
% and column of C, so one extra coordinate and a duplicated last row/column
% are added. This shifts each cell boundary by half a step, which is not
% visible at the cell counts plotted here. NaN cells are left undrawn.
%
% INPUTS
%   ax  target axes
%   x   [1 x nCol] or [nCol x 1], column (e.g. time) coordinates
%   y   [1 x nRow] or [nRow x 1], row (e.g. depth) coordinates
%   C   [nRow x nCol], same orientation imagesc(ax,x,y,C) expects
%
% OUTPUT  h, the pcolor surface handle.

maxCells = 1200000;
assert(numel(C) <= maxCells, 'pcolorLikeImagesc:tooManyCells', ...
    'C has %d cells, over the %d-cell limit; rebin the data for display first.', numel(C), maxCells);

x = x(:)'; y = y(:)';
dx = median(diff(x), 'omitnan'); if ~isfinite(dx) || dx == 0; dx = 1; end
dy = median(diff(y), 'omitnan'); if ~isfinite(dy) || dy == 0; dy = 1; end
xEdges = [x, x(end)+dx];
yEdges = [y, y(end)+dy];
Cpad = [C, C(:,end); C(end,:), C(end,end)];

h = pcolor(ax, xEdges, yEdges, Cpad);
set(h, 'EdgeColor', 'none');
end
