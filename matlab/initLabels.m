function labels = initLabels(windows)
%INITLABELS Empty label struct array, one element per window (see
% buildLabelWindows.m), for labelMSLTopEdgeGUI.m.
%
% OUTPUT  labels, struct array with fields:
%   deploymentID, winStart, winEnd  echoed from windows
%   pointsTime       [n x 1] datetime, hand-clicked points along the true top
%                    edge, ascending in time (empty until labeled). May have
%                    gaps in time where the layer disappears and reappears
%                    within the window.
%   pointsDepth      [n x 1] double, matching depths (m)
%   isAbsent         logical, true if the whole window was marked "no clear
%                    boundary anywhere"
%   absentIntervals  [m x 2] datetime, each row [start, end] a sub-range of
%                    this window with no clear boundary, alongside points
%                    traced elsewhere in the same window. Time covered by
%                    neither points nor an absent interval is unlabeled and
%                    skipped by the scoring.
%   labeledAt        datetime this window was last saved, NaT if untouched

n = numel(windows);
labels = struct('deploymentID', {windows.deploymentID}, 'winStart', {windows.winStart}, ...
    'winEnd', {windows.winEnd}, 'pointsTime', cell(1,n), 'pointsDepth', cell(1,n), ...
    'isAbsent', num2cell(false(1,n)), 'absentIntervals', cell(1,n), 'labeledAt', cell(1,n));
for i = 1:n
    labels(i).pointsTime = datetime.empty(0,1);
    labels(i).pointsDepth = zeros(0,1);
    labels(i).absentIntervals = datetime.empty(0,2);
    labels(i).labeledAt = NaT;
end
end
