function edge = computeMSLEdgeAndCore(nasc, cfg)
%COMPUTEMSLEDGEANDCORE Top edge and core of the non-migratory mesopelagic
% scattering layer (MSL) in an echosounder NASC depth x time grid.
%
% The detector seeds on the layer core and walks toward the surface with
% hysteresis. For every time step:
%   1. Reference: the quiet water in cfg.edgeRefDepthRangeM (shallower than
%      the search range). Robust median and MAD give z = (NASC - median) / std.
%   2. Core: the maximum of the smoothed profile within
%      cfg.edgeSearchDepthRangeM. It counts as a layer only if it clears
%      z >= cfg.edgeZHigh and an absolute floor (the
%      cfg.edgeMinAbsNASCPercentile-th percentile of positive in-window NASC,
%      pooled over the record). Otherwise the step has no layer.
%   3. Walk: from the core, step toward the surface one cell at a time while
%      z >= cfg.edgeZLow and NASC >= floor * cfg.edgeAbsFloorLowFrac. The top
%      edge is the sub-bin interpolated crossing of edgeZLow where the walk stops.
%   4. Core output: depth and raw NASC of the seed cell. The edge-to-core gap
%      is not constrained.
%   5. Integrated NASC: raw NASC summed from the top edge to the bottom of
%      the search range, times the depth step.
%
% Time tracking: a continuity gate around the last edge (base width plus a
% widening rate per hour), a timeout after which the track is declared lost,
% and a consecutive-detection requirement to re-acquire a lost track.
%
% INPUTS
%   nasc  struct: time [nTime x 1 datetime], depth [1 x nDepth] metres below
%         the surface, ASCENDING, NASC [nTime x nDepth] linear NASC (NaN is
%         treated as 0). Optional realTime, nightIndex and deploymentID
%         [nTime x 1] pass through; a change in deploymentID resets tracking.
%   cfg   struct with edgeSearchDepthRangeM [min max], edgeRefDepthRangeM
%         [min max], edgeZHigh, edgeZLow, edgeMinAbsNASCPercentile and
%         edgeAbsFloorLowFrac. Optional, with defaults: nascLayerMinThicknessM (20),
%         edgeContinuityGateBaseM (20), edgeContinuityGateRateMPerHour (4),
%         edgeAbsenceMinDurationHours (8), edgeReacquireMinConsecutive (1),
%         edgeReacquireToleranceM (= edgeContinuityGateBaseM).
%
% OUTPUT  struct edge, all [nTime x 1] aligned to nasc.time: edgeDepthTop,
%   edgeDepthTopSharpness, coreDepth, coreAmplitude, integratedNASC,
%   trackState, layerPresent, refMedian, refRobustStd, plus time, realTime,
%   nightIndex, deploymentID and params.
depthVec = nasc.depth(:)';
assert(issorted(depthVec), 'computeMSLEdgeAndCore:depthNotSorted', 'nasc.depth must be ascending.');
depthStep = median(diff(depthVec));

searchMin = cfg.edgeSearchDepthRangeM(1);
searchMax = cfg.edgeSearchDepthRangeM(2);
refMin = cfg.edgeRefDepthRangeM(1);
refMax = cfg.edgeRefDepthRangeM(2);
assert(refMax <= searchMin, 'computeMSLEdgeAndCore:refNotShallower', ...
    'cfg.edgeRefDepthRangeM must be shallower than (or touch) cfg.edgeSearchDepthRangeM(1).');

searchIdx = find(depthVec >= searchMin & depthVec <= searchMax);
refIdx = find(depthVec >= refMin & depthVec <= refMax);
assert(~isempty(searchIdx), 'computeMSLEdgeAndCore:emptySearch', 'Search range has no overlap with nasc.depth.');
assert(~isempty(refIdx), 'computeMSLEdgeAndCore:emptyRef', 'Reference range has no overlap with nasc.depth.');

minThicknessM = getOr(cfg, 'nascLayerMinThicknessM', 20);
smoothWindowCells = max(1, round((minThicknessM/2) / depthStep));

zHigh = getOr(cfg, 'edgeZHigh', 4);
zLow = getOr(cfg, 'edgeZLow', 2);
minAbsPctile = getOr(cfg, 'edgeMinAbsNASCPercentile', 50);
absFloorLowFrac = getOr(cfg, 'edgeAbsFloorLowFrac', 0.5);

gateBaseM = getOr(cfg, 'edgeContinuityGateBaseM', 20);
gateRateMPerHour = getOr(cfg, 'edgeContinuityGateRateMPerHour', 4);
absenceMinDurationHours = getOr(cfg, 'edgeAbsenceMinDurationHours', 8);
reacquireMinConsecutive = getOr(cfg, 'edgeReacquireMinConsecutive', 1);
reacquireToleranceM = getOr(cfg, 'edgeReacquireToleranceM', gateBaseM);

nascFilled = nasc.NASC;
nascFilled(~isfinite(nascFilled)) = 0;
nascSmooth = movmean(nascFilled, smoothWindowCells, 2);

% Absolute-magnitude floor, pooled over the whole record: a near-zero
% reference would make a pure z-test degenerate.
nascForFloor = nascSmooth(:, searchIdx);
nascForFloor = nascForFloor(isfinite(nascForFloor) & nascForFloor > 0);
minAbsNASC = prctile(nascForFloor, minAbsPctile);
fprintf('computeMSLEdgeAndCore: absolute floor = %.4f (%g-th percentile of positive in-search-window NASC).\n', ...
    minAbsNASC, minAbsPctile);

if isfield(nasc, 'deploymentID') && ~isempty(nasc.deploymentID)
    deployID = nasc.deploymentID(:);
else
    deployID = repmat("ALL", numel(nasc.time), 1);
end
segmentBreak = [true; deployID(2:end) ~= deployID(1:end-1)];
timeHours = hours(nasc.time(:) - nasc.time(1));

nTime = numel(nasc.time);
edgeDepthTop = nan(nTime,1); edgeDepthTopSharpness = nan(nTime,1);
coreDepthOut = nan(nTime,1); coreAmplitudeOut = nan(nTime,1); integratedNASCOut = nan(nTime,1);
refMedianOut = nan(nTime,1); refRobustStdOut = nan(nTime,1);
trackStateStr = repmat("lost", nTime, 1);

state = 'lost'; anchorDepth = nan; anchorTimeHours = nan; failHours = 0;
pendingDepth = nan; pendingCount = 0;

for it = 1:nTime
    if segmentBreak(it)
        state = 'lost'; anchorDepth = nan; anchorTimeHours = nan; failHours = 0; pendingDepth = nan; pendingCount = 0;
    end
    dtHours = 0;
    if it > 1 && ~segmentBreak(it)
        dtHours = max(timeHours(it) - timeHours(it-1), 0);
    end

    row = nascSmooth(it, :);
    refVals = row(refIdx);
    refMedian = median(refVals);
    refRobustStd = max(1.4826 * mad(refVals, 1), 1e-6);
    refMedianOut(it) = refMedian; refRobustStdOut(it) = refRobustStd;

    zRow = (row - refMedian) / refRobustStd;

    found = false; crossDepth = nan; sharp = nan; coreDepth = nan; coreAmplitude = nan; integratedNASC = nan;

    % SEED: the profile's own global max within the search range --
    % unambiguous regardless of how quiet or busy the water above it is.
    [coreValSmooth, coreRelIdx] = max(row(searchIdx));
    coreGlobalIdx = searchIdx(coreRelIdx);

    if zRow(coreGlobalIdx) >= zHigh && coreValSmooth >= minAbsNASC
        % WALK toward shallow from the seed through the LOW (hysteresis)
        % tier, stopping at the first cell that fails it.
        minAbsNASCLow = minAbsNASC * absFloorLowFrac;
        walkIdx = coreGlobalIdx;
        for gi = (coreGlobalIdx-1):-1:searchIdx(1)
            if zRow(gi) >= zLow && row(gi) >= minAbsNASCLow
                walkIdx = gi;
            else
                break
            end
        end

        prevGlobalIdx = walkIdx - 1;
        if prevGlobalIdx >= 1 && ~(zRow(prevGlobalIdx) >= zLow && row(prevGlobalIdx) >= minAbsNASCLow)
            z0 = zRow(prevGlobalIdx); z1 = zRow(walkIdx);
            d0 = depthVec(prevGlobalIdx); d1 = depthVec(walkIdx);
            frac = min(max((zLow - z0) / max(z1 - z0, eps), 0), 1);
            crossDepth = d0 + frac*(d1-d0);
            sharp = abs(z1-z0) / max(d1-d0, eps);
        else
            crossDepth = depthVec(walkIdx);
        end
        found = true;

        % CORE: the seed's own depth and RAW (unsmoothed) NASC value.
        coreDepth = depthVec(coreGlobalIdx);
        coreAmplitude = nascFilled(it, coreGlobalIdx);
        coreSearchIdx = searchIdx(depthVec(searchIdx) >= crossDepth);
        integratedNASC = sum(nascFilled(it, coreSearchIdx)) * depthStep;
    end

    wasLost = strcmp(state, 'lost');
    if found && ~wasLost
        dtHoursSinceAnchor = timeHours(it) - anchorTimeHours;
        gateM = min(searchMax - searchMin, gateBaseM + gateRateMPerHour * dtHoursSinceAnchor);
        if abs(crossDepth - anchorDepth) > gateM
            found = false;
        end
    end

    if found && wasLost
        if pendingCount > 0 && abs(crossDepth - pendingDepth) <= reacquireToleranceM
            pendingCount = pendingCount + 1;
        else
            pendingCount = 1;
        end
        pendingDepth = crossDepth;
        if pendingCount < reacquireMinConsecutive
            found = false;
        else
            pendingDepth = nan; pendingCount = 0;
        end
    elseif ~found
        pendingDepth = nan; pendingCount = 0;
    end

    if found
        edgeDepthTop(it) = crossDepth;
        edgeDepthTopSharpness(it) = sharp;
        coreDepthOut(it) = coreDepth; coreAmplitudeOut(it) = coreAmplitude;
        integratedNASCOut(it) = integratedNASC;
        anchorDepth = crossDepth; anchorTimeHours = timeHours(it);
        state = 'tracking'; trackStateStr(it) = "found"; failHours = 0;
    else
        trackStateStr(it) = "lost";
        if strcmp(state, 'tracking')
            failHours = failHours + max(dtHours, 0);
            if failHours > absenceMinDurationHours
                state = 'lost'; anchorDepth = nan; anchorTimeHours = nan;
            end
        end
    end
end

edge.time = nasc.time;
if isfield(nasc, 'realTime') && ~isempty(nasc.realTime)
    edge.realTime = nasc.realTime;
end
if isfield(nasc, 'nightIndex') && ~isempty(nasc.nightIndex)
    edge.nightIndex = nasc.nightIndex;
end
edge.deploymentID = deployID;
edge.edgeDepthTop = edgeDepthTop;
edge.edgeDepthTopSharpness = edgeDepthTopSharpness;
edge.coreDepth = coreDepthOut;
edge.coreAmplitude = coreAmplitudeOut;
edge.integratedNASC = integratedNASCOut;
edge.refMedian = refMedianOut; edge.refRobustStd = refRobustStdOut;
edge.trackState = categorical(trackStateStr);
edge.layerPresent = isfinite(edgeDepthTop);
edge.params = struct('searchDepthRangeM', [searchMin searchMax], 'refDepthRangeM', [refMin refMax], ...
    'zHigh', zHigh, 'zLow', zLow, 'minAbsNASCPercentile', minAbsPctile, 'absFloorLowFrac', absFloorLowFrac, ...
    'gateBaseM', gateBaseM, 'gateRateMPerHour', gateRateMPerHour, 'absenceMinDurationHours', absenceMinDurationHours);

nPresent = sum(edge.layerPresent);
fprintf('computeMSLEdgeAndCore: %d/%d steps tracked (%.1f%%). Edge %.0f-%.0f m, core %.0f-%.0f m when tracked.\n', ...
    nPresent, nTime, 100*nPresent/nTime, min(edgeDepthTop,[],'omitnan'), max(edgeDepthTop,[],'omitnan'), ...
    min(coreDepthOut,[],'omitnan'), max(coreDepthOut,[],'omitnan'));
end

function v = getOr(cfg, field, default)
if isfield(cfg, field) && ~isempty(cfg.(field))
    v = cfg.(field);
else
    v = default;
end
end
