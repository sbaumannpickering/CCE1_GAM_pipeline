function windows = buildLabelWindows(nascNative, cfg)
%BUILDLABELWINDOWS Evenly spaced time windows for hand-labeling the layer top
% edge. cfg.edgeLabelWindowsPerDeployment windows of cfg.edgeLabelWindowHours
% each, spaced evenly across each deployment's span. Windows are chosen by time
% only, with no reference to the detector output, so the label set is
% independent of the detector it is used to score.
%
% INPUTS
%   nascNative  NASC grid struct (time, depth, NASC, optional deploymentID)
%   cfg         struct with edgeLabelWindowHours (default 48) and
%               edgeLabelWindowsPerDeployment (default 8)
%
% OUTPUT  struct array with fields deploymentID, winStart, winEnd, ordered by
%         deployment then start time.

nascNative = ensureDeploymentID(nascNative);
windowHours = getOr(cfg, 'edgeLabelWindowHours', 48);
nPerDeploy = getOr(cfg, 'edgeLabelWindowsPerDeployment', 8);

uDeploy = unique(nascNative.deploymentID(nascNative.deploymentID ~= ""), 'stable');
windows = struct('deploymentID', {}, 'winStart', {}, 'winEnd', {});
for d = 1:numel(uDeploy)
    mask = nascNative.deploymentID == uDeploy(d);
    tDep = nascNative.time(mask);
    t0 = tDep(1); t1 = tDep(end);
    usableSpanHours = hours(t1 - t0) - windowHours;
    assert(usableSpanHours > 0, 'buildLabelWindows:deploymentTooShort', ...
        '%s spans only %.1f h, shorter than one %g-h label window.', uDeploy(d), hours(t1-t0), windowHours);
    if nPerDeploy == 1
        startsHours = 0;
    else
        startsHours = linspace(0, usableSpanHours, nPerDeploy);
    end
    for w = 1:numel(startsHours)
        ws = t0 + hours(startsHours(w));
        windows(end+1).deploymentID = uDeploy(d); %#ok<AGROW>
        windows(end).winStart = ws;
        windows(end).winEnd = ws + hours(windowHours);
    end
end

fprintf('buildLabelWindows: %d window(s) (%g h each) across %d deployment(s), %d per deployment.\n', ...
    numel(windows), windowHours, numel(uDeploy), nPerDeploy);
end

function v = getOr(cfg, field, default)
if isfield(cfg, field) && ~isempty(cfg.(field))
    v = cfg.(field);
else
    v = default;
end
end
