function labelMSLTopEdgeGUI(nascStitched, cfg, windowsOverride, labelFileOverride)
%LABELMSLTOPEDGEGUI Interactive tool to hand-label the top edge of the
% mesopelagic scattering layer (MSL) on the stitched night-time grid. The
% labels are the ground truth that scoreMSLEdgeAgainstLabels.m compares the
% detector against.
%
% The grid comes from restrictToNightTrimmed.m: nights are joined back to back
% on a synthetic, gap-free time axis. Everything clicked is stored in that
% stitched time. The x-axis tick labels, the title and the grey dotted vertical
% lines show the real calendar date and the real night boundaries.
%
% CONTROLS
%   TRACE mode (default)  left-click adds a point on the true top edge (points
%                         in a window are connected); right-click removes the
%                         nearest point.
%   MARK ABSENT mode      (toggle button) left-click-drag shades a span with no
%                         clear boundary; right-click removes the nearest band.
%   << Prev / Next >>     navigate windows (saves on every move)
%   Undo point            remove the most recent click
%   Clear window          remove everything in the current window
%   Mark ABSENT / no clear layer  mark the whole window as having no layer, then
%                         advance
%   Save now              write progress immediately (also on navigating and on
%                         closing the window)
% Running the tool again resumes at the first window that is not yet labeled.
% If the window plan changes, labels are carried over by matching window start
% time and the old file is backed up first.
%
% INPUTS
%   nascStitched  stitched grid from restrictToNightTrimmed.m, with time, depth,
%                 NASC, realTime and nightIndex
%   cfg           struct with edgeSearchDepthRangeM, edgeLabelWindowHours,
%                 edgeLabelWindowsPerDeployment and edgeLabelFile
%   windowsOverride, labelFileOverride  optional: caller-supplied windows
%                 (struct array as from buildLabelWindows.m) and label file path
%
% Needs an interactive MATLAB session with a display (not matlab -batch).

assert(isfield(nascStitched, 'realTime') && isfield(nascStitched, 'nightIndex'), ...
    'labelMSLTopEdgeGUI:missingFields', ...
    'nascStitched must have .realTime and .nightIndex (from restrictToNightTrimmed.m).');
nascStitched = ensureDeploymentID(nascStitched);

if nargin < 3 || isempty(windowsOverride)
    windowsOverride = [];
end
if nargin < 4 || isempty(labelFileOverride)
    labelFileOverride = cfg.edgeLabelFile;
end
labelFile = labelFileOverride;

if isempty(windowsOverride)
    windows = buildLabelWindows(nascStitched, cfg);
else
    windows = windowsOverride(:)';
    fprintf('labelMSLTopEdgeGUI: using %d caller-supplied window(s); labels -> %s\n', numel(windows), labelFile);
end

if isfile(labelFile)
    S = load(labelFile, 'labels');
    labels = S.labels;
    if numel(labels) ~= numel(windows)
        backupFile = [labelFile(1:end-4) '_backup_' char(datetime('now', 'Format', 'yyyyMMdd_HHmmss')) '.mat'];
        copyfile(labelFile, backupFile);
        if ~isfield(labels, 'absentIntervals')
            for i = 1:numel(labels)
                labels(i).absentIntervals = datetime.empty(0,2);
            end
        end
        oldLabels = labels;
        labels = initLabels(windows);
        winKey = @(s) string(s.deploymentID) + "|" + string(s.winStart, 'yyyyMMddHHmmss');
        oldKeys = arrayfun(winKey, oldLabels);
        nCarried = 0;
        for i = 1:numel(labels)
            j = find(oldKeys == winKey(labels(i)), 1);
            if ~isempty(j)
                for f = string(fieldnames(labels))'
                    labels(i).(f) = oldLabels(j).(f);
                end
                nCarried = nCarried + 1;
            end
        end
        warning('labelMSLTopEdgeGUI:windowCountChanged', ...
            ['Saved label file has %d window(s) but the current window plan has %d -- carried over %d matching ' ...
             'label(s), %d of the saved ones no longer match a window and were dropped. Backup of the old file: %s'], ...
            numel(oldLabels), numel(windows), nCarried, numel(oldLabels) - nCarried, backupFile);
    else
        if ~isfield(labels, 'absentIntervals')
            for i = 1:numel(labels)
                labels(i).absentIntervals = datetime.empty(0,2);
            end
        end
        fprintf('Resuming from %s (%d/%d window(s) already labeled).\n', labelFile, ...
            sum(~cellfun(@isempty, {labels.pointsTime}) | [labels.isAbsent]), numel(labels));
    end
else
    labels = initLabels(windows);
end

startIdx = find(cellfun(@isempty, {labels.pointsTime}) & ~[labels.isAbsent], 1);
if isempty(startIdx)
    startIdx = 1;
end

fig = figure('Name', 'MSL Top-Edge Labeling (stitched nights, real dates on axis)', 'NumberTitle', 'off', ...
    'Units', 'pixels', 'Position', [50 50 1700 920], 'Color', 'w', 'CloseRequestFcn', @onClose, ...
    'MenuBar', 'none', 'ToolBar', 'figure');
ax = axes(fig, 'Units', 'normalized', 'Position', [0.055 0.17 0.92 0.76]);

uicontrol(fig, 'Style', 'text', 'Units', 'pixels', 'Position', [20 882 1660 30], 'HorizontalAlignment', 'left', ...
    'FontSize', 10, 'BackgroundColor', 'w', 'String', ...
    ['TRACE mode: left-click adds a boundary point (= where YOU think the MSL top edge is), right-click removes nearest. ' ...
     'MARK-ABSENT mode: left-click-drag selects a no-clear-boundary span. Grey dotted lines = real night boundaries.']);

btnY = 52; btnH = 32;
uicontrol(fig, 'Style', 'pushbutton', 'Units', 'pixels', 'Position', [20 btnY 90 btnH], 'String', '<< Prev', 'Callback', @(s,e) navigate(fig, -1));
uicontrol(fig, 'Style', 'pushbutton', 'Units', 'pixels', 'Position', [120 btnY 90 btnH], 'String', 'Next >>', 'Callback', @(s,e) navigate(fig, 1));
uicontrol(fig, 'Style', 'pushbutton', 'Units', 'pixels', 'Position', [230 btnY 110 btnH], 'String', 'Undo point', 'Callback', @(s,e) undoPoint(fig));
uicontrol(fig, 'Style', 'pushbutton', 'Units', 'pixels', 'Position', [350 btnY 120 btnH], 'String', 'Clear window', 'Callback', @(s,e) clearWindow(fig));
uicontrol(fig, 'Style', 'pushbutton', 'Units', 'pixels', 'Position', [480 btnY 220 btnH], 'String', 'Mark ABSENT / no clear layer', 'Callback', @(s,e) markAbsent(fig));
modeBtn = uicontrol(fig, 'Style', 'pushbutton', 'Units', 'pixels', 'Position', [710 btnY 300 btnH], 'Callback', @(s,e) toggleMode(fig));
uicontrol(fig, 'Style', 'pushbutton', 'Units', 'pixels', 'Position', [1020 btnY 100 btnH], 'String', 'Save now', 'Callback', @(s,e) saveLabels(fig, true));
progressText = uicontrol(fig, 'Style', 'text', 'Units', 'pixels', 'Position', [20 12 1660 34], ...
    'HorizontalAlignment', 'left', 'FontSize', 10, 'BackgroundColor', 'w');

state.nascNative = nascStitched;
state.windows = windows;
state.labels = labels;
state.idx = startIdx;
state.ax = ax;
state.progressText = progressText;
state.modeBtn = modeBtn;
state.cfg = cfg;
state.labelFile = labelFile;
state.dirty = false;
state.mode = 'trace';
state.dragging = false;
state.dragStartT = NaT;
state.dragPatch = gobjects(0);
guidata(fig, state);

fig.Units = 'pixels';
fig.WindowButtonDownFcn = @(s, e) onClick(fig);
fig.WindowButtonMotionFcn = @(s, e) onMotion(fig);
fig.WindowButtonUpFcn = @(s, e) onRelease(fig);

updateModeButton(fig);
renderWindow(fig);
end

%% ---- callbacks / helpers ----

function renderWindow(fig)
state = guidata(fig);
w = state.windows(state.idx);
nascNative = state.nascNative;
ax = state.ax;
cfg = state.cfg;

cla(ax); hold(ax, 'on');
depthPad = [cfg.edgeSearchDepthRangeM(1)-50, cfg.edgeSearchDepthRangeM(2)+50];
inCrop = nascNative.depth >= depthPad(1) & nascNative.depth <= depthPad(2);
depthCropped = nascNative.depth(inCrop);
inWin = nascNative.time >= w.winStart & nascNative.time <= w.winEnd & nascNative.deploymentID == w.deploymentID;
tWin = nascNative.time(inWin);
nascWin = nascNative.NASC(inWin, inCrop);

renderEchogramLabelStyle(ax, tWin, depthCropped, nascWin);
yline(ax, cfg.edgeSearchDepthRangeM(1), ':', 'Color', [0.85 0.85 0.85], 'HandleVisibility', 'off');
yline(ax, cfg.edgeSearchDepthRangeM(2), ':', 'Color', [0.85 0.85 0.85], 'HandleVisibility', 'off');

xlim(ax, [w.winStart, w.winEnd]);
ylim(ax, depthPad);

% --- real-date tick labels and night-boundary lines (nascNative.time is a
% stitched pseudo-time axis, so its own values are not real dates) ---
realTimeWin = nascNative.realTime(inWin);
nightIdxWin = nascNative.nightIndex(inWin);
nTicks = min(8, numel(tWin));
if nTicks >= 2
    tickPos = round(linspace(1, numel(tWin), nTicks));
    ax.XTick = tWin(tickPos);
    ax.XTickLabel = cellstr(datestr(realTimeWin(tickPos), 'dd-mmm HH:MM')); %#ok<DATST>
    xtickangle(ax, 30);
end
nightBoundaryPos = find([true; diff(nightIdxWin) ~= 0]);
for b = nightBoundaryPos(:)'
    xline(ax, tWin(b), ':', 'Color', [0.6 0.6 0.6], 'HandleVisibility', 'off');
end

ylabel(ax, 'Depth (m)');
xlabel(ax, 'Time (stitched -- no day gaps; tick labels & grey lines show REAL date/night boundaries)');
if ~isempty(realTimeWin)
    realRangeStr = sprintf('REAL: %s to %s', datestr(realTimeWin(1), 'dd-mmm-yyyy HH:MM'), ...
        datestr(realTimeWin(end), 'dd-mmm-yyyy HH:MM')); %#ok<DATST>
else
    realRangeStr = '';
end
title(ax, sprintf('[%d/%d]  %s   %s', state.idx, numel(state.windows), w.deploymentID, realRangeStr), ...
    'Interpreter', 'none');

lbl = state.labels(state.idx);
if lbl.isAbsent
    text(ax, w.winStart + (w.winEnd-w.winStart)/2, mean(depthPad), 'MARKED ABSENT -- no clear boundary anywhere', ...
        'Color', [1 0.2 0.2], 'FontSize', 18, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
else
    for k = 1:size(lbl.absentIntervals, 1)
        patch(ax, [lbl.absentIntervals(k,1), lbl.absentIntervals(k,1), lbl.absentIntervals(k,2), lbl.absentIntervals(k,2)], ...
            [depthPad(1), depthPad(2), depthPad(2), depthPad(1)], [0.65 0.1 0.1], ...
            'FaceAlpha', 0.30, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    end
    if ~isempty(lbl.pointsTime)
        plot(ax, lbl.pointsTime, lbl.pointsDepth, '-', 'Color', [1 1 1], 'LineWidth', 2.5, 'HandleVisibility', 'off');
        plot(ax, lbl.pointsTime, lbl.pointsDepth, 'o-', 'Color', [1 0.1 0.1], 'MarkerFaceColor', [1 1 0.2], ...
            'MarkerEdgeColor', [0 0 0], 'LineWidth', 1.5, 'MarkerSize', 7, 'HandleVisibility', 'off');
    end
end

updateProgressText(fig);
end

function updateProgressText(fig)
state = guidata(fig);
nAbsentWhole = sum([state.labels.isAbsent]);
nPoints = sum(~cellfun(@isempty, {state.labels.pointsTime}));
nPartial = sum(cellfun(@(x) size(x,1)>0, {state.labels.absentIntervals}));
nDone = sum(~cellfun(@isempty, {state.labels.pointsTime}) | [state.labels.isAbsent]);
dirtyFlag = '';
if state.dirty
    dirtyFlag = '  [unsaved changes]';
end
state.progressText.String = sprintf(['Window %d of %d.  %d labeled (%d traced [%d with a partial absent band], %d marked fully absent), ' ...
    '%d remaining.%s\nFile: %s'], state.idx, numel(state.windows), nDone, nPoints, nPartial, nAbsentWhole, ...
    numel(state.windows)-nDone, dirtyFlag, state.labelFile);
end

function toggleMode(fig)
state = guidata(fig);
if strcmp(state.mode, 'trace')
    state.mode = 'absent';
else
    state.mode = 'trace';
end
guidata(fig, state);
updateModeButton(fig);
end

function updateModeButton(fig)
state = guidata(fig);
if strcmp(state.mode, 'trace')
    state.modeBtn.String = 'Mode: TRACE boundary (click to switch)';
    state.modeBtn.BackgroundColor = [0.85 0.93 1];
else
    state.modeBtn.String = 'Mode: MARK ABSENT range (click to switch)';
    state.modeBtn.BackgroundColor = [1 0.8 0.8];
end
end

function [tClick, dClick, inAxes] = clickToData(fig)
state = guidata(fig);
ax = state.ax;
axPos = getpixelposition(ax, true);
cp = fig.CurrentPoint;
fracX = (cp(1) - axPos(1)) / axPos(3);
fracY = (cp(2) - axPos(2)) / axPos(4);
inAxes = fracX >= 0 && fracX <= 1 && fracY >= 0 && fracY <= 1;
fracXClamped = min(max(fracX, 0), 1);
fracYClamped = min(max(fracY, 0), 1);
xl = xlim(ax); yl = ylim(ax);
tClick = xl(1) + fracXClamped * (xl(2) - xl(1));
dClick = yl(2) - fracYClamped * (yl(2) - yl(1));
end

function onClick(fig)
[tClick, dClick, inAxes] = clickToData(fig);
if ~inAxes
    return
end
state = guidata(fig);

if strcmp(state.mode, 'trace')
    switch fig.SelectionType
        case 'normal'
            addPoint(fig, tClick, dClick);
        case 'alt'
            removeNearestPoint(fig, tClick, dClick);
    end
else
    switch fig.SelectionType
        case 'normal'
            ax = state.ax;
            yl = ylim(ax);
            state.dragging = true;
            state.dragStartT = tClick;
            state.dragPatch = patch(ax, [tClick tClick tClick tClick], [yl(1) yl(2) yl(2) yl(1)], [0.65 0.1 0.1], ...
                'FaceAlpha', 0.30, 'EdgeColor', [0.4 0 0], 'LineWidth', 1, 'HandleVisibility', 'off');
            guidata(fig, state);
        case 'alt'
            removeNearestAbsentInterval(fig, tClick);
    end
end
end

function onMotion(fig)
state = guidata(fig);
if ~state.dragging
    return
end
[tNow, ~, ~] = clickToData(fig);
tLo = min(state.dragStartT, tNow); tHi = max(state.dragStartT, tNow);
if isvalid(state.dragPatch)
    state.dragPatch.XData = [tLo, tLo, tHi, tHi];
end
end

function onRelease(fig)
state = guidata(fig);
if ~state.dragging
    return
end
[tNow, ~, ~] = clickToData(fig);
tLo = min(state.dragStartT, tNow); tHi = max(state.dragStartT, tNow);
state.dragging = false;
if isvalid(state.dragPatch)
    delete(state.dragPatch);
end
guidata(fig, state);
if hours(tHi - tLo) > 0.02
    addAbsentInterval(fig, tLo, tHi);
else
    renderWindow(fig);
end
end

function addPoint(fig, tClick, dClick)
state = guidata(fig);
lbl = state.labels(state.idx);
lbl.isAbsent = false;
lbl.pointsTime(end+1, 1) = tClick;
lbl.pointsDepth(end+1, 1) = dClick;
[lbl.pointsTime, order] = sort(lbl.pointsTime);
lbl.pointsDepth = lbl.pointsDepth(order);
state.labels(state.idx) = lbl;
state.dirty = true;
guidata(fig, state);
renderWindow(fig);
end

function removeNearestPoint(fig, tClick, dClick)
state = guidata(fig);
lbl = state.labels(state.idx);
if isempty(lbl.pointsTime)
    return
end
ax = state.ax;
xl = xlim(ax); yl = ylim(ax);
tScale = hours(xl(2) - xl(1)); dScale = diff(yl);
dt = hours(lbl.pointsTime - tClick) / tScale;
dd = (lbl.pointsDepth - dClick) / dScale;
[~, nearest] = min(dt.^2 + dd.^2);
lbl.pointsTime(nearest) = [];
lbl.pointsDepth(nearest) = [];
state.labels(state.idx) = lbl;
state.dirty = true;
guidata(fig, state);
renderWindow(fig);
end

function addAbsentInterval(fig, tLo, tHi)
state = guidata(fig);
lbl = state.labels(state.idx);
lbl.isAbsent = false;
lbl.absentIntervals = mergeIntervals([lbl.absentIntervals; tLo, tHi]);
state.labels(state.idx) = lbl;
state.dirty = true;
guidata(fig, state);
renderWindow(fig);
end

function removeNearestAbsentInterval(fig, tClick)
state = guidata(fig);
lbl = state.labels(state.idx);
if isempty(lbl.absentIntervals)
    return
end
starts = lbl.absentIntervals(:,1); ends = lbl.absentIntervals(:,2);
inside = tClick >= starts & tClick <= ends;
if any(inside)
    remove = find(inside, 1);
else
    distToEdge = min(abs(hours(tClick - starts)), abs(hours(tClick - ends)));
    [~, remove] = min(distToEdge);
end
lbl.absentIntervals(remove, :) = [];
state.labels(state.idx) = lbl;
state.dirty = true;
guidata(fig, state);
renderWindow(fig);
end

function merged = mergeIntervals(intervals)
if isempty(intervals)
    merged = intervals;
    return
end
[~, order] = sort(intervals(:,1));
intervals = intervals(order, :);
merged = intervals(1,:);
for i = 2:size(intervals,1)
    if intervals(i,1) <= merged(end,2)
        merged(end,2) = max(merged(end,2), intervals(i,2));
    else
        merged(end+1,:) = intervals(i,:); %#ok<AGROW>
    end
end
end

function undoPoint(fig)
state = guidata(fig);
lbl = state.labels(state.idx);
if ~isempty(lbl.pointsTime)
    lbl.pointsTime(end) = [];
    lbl.pointsDepth(end) = [];
    state.labels(state.idx) = lbl;
    state.dirty = true;
    guidata(fig, state);
    renderWindow(fig);
end
end

function clearWindow(fig)
state = guidata(fig);
lbl = state.labels(state.idx);
lbl.pointsTime = datetime.empty(0,1);
lbl.pointsDepth = zeros(0,1);
lbl.absentIntervals = datetime.empty(0,2);
lbl.isAbsent = false;
state.labels(state.idx) = lbl;
state.dirty = true;
guidata(fig, state);
renderWindow(fig);
end

function markAbsent(fig)
state = guidata(fig);
lbl = state.labels(state.idx);
lbl.pointsTime = datetime.empty(0,1);
lbl.pointsDepth = zeros(0,1);
lbl.absentIntervals = datetime.empty(0,2);
lbl.isAbsent = true;
lbl.labeledAt = datetime('now');
state.labels(state.idx) = lbl;
state.dirty = true;
guidata(fig, state);
saveLabels(fig, false);
navigate(fig, 1);
end

function navigate(fig, delta)
state = guidata(fig);
lbl = state.labels(state.idx);
if ~isempty(lbl.pointsTime) || lbl.isAbsent || ~isempty(lbl.absentIntervals)
    lbl.labeledAt = datetime('now');
    state.labels(state.idx) = lbl;
    guidata(fig, state);
end
saveLabels(fig, false);
state = guidata(fig);
newIdx = state.idx + delta;
if newIdx < 1 || newIdx > numel(state.windows)
    return
end
state.idx = newIdx;
guidata(fig, state);
renderWindow(fig);
end

function saveLabels(fig, verbose)
state = guidata(fig);
labels = state.labels; %#ok<NASGU>
save(state.labelFile, 'labels');
state.dirty = false;
guidata(fig, state);
if verbose
    fprintf('Saved %d window label(s) to %s\n', numel(state.labels), state.labelFile);
end
if isvalid(state.progressText)
    updateProgressText(fig);
end
end

function onClose(fig, ~)
saveLabels(fig, true);
delete(fig);
end
