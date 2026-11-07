function score = scoreMSLEdgeAgainstLabels(nascStitched, cfg, labelFileOverride, detectorFcn)
%SCOREMSLEDGEAGAINSTLABELS Error of a top-edge detector against the hand
% labels made with labelMSLTopEdgeGUI.m.
%
% The detector is re-run on nascStitched with cfg, so changing cfg and scoring
% again tests a parameter change. Labels and detector output are both in the
% stitched pseudo-time of restrictToNightTrimmed.m. Any detector that takes
% (nascStitched, cfg) and returns a struct with .time, .deploymentID and
% .edgeDepthTop works.
%
% Traced windows give a depth error and a miss rate (time inside "absent"
% bands is excluded). Absent intervals give a false-positive rate (detector
% steps that report a layer where the labeler marked none).
%
% INPUTS
%   nascStitched       the stitched grid the labels were traced on
%   cfg                struct with the detector's fields and edgeLabelFile
%   labelFileOverride  optional label file (default cfg.edgeLabelFile)
%   detectorFcn        optional function handle (default @computeMSLEdgeAndCore)
%
% OUTPUT  struct score with perWindow, absentIntervals, overallMeanAbsErrorM,
%   overallBiasM (positive = detector too deep), overallRmsErrorM,
%   overallMissRate, overallFalsePositiveRate, nTracedWindows and
%   nAbsentIntervals.

nascStitched = ensureDeploymentID(nascStitched);
if nargin < 3 || isempty(labelFileOverride)
    labelFile = cfg.edgeLabelFile;
else
    labelFile = labelFileOverride;
end
if nargin < 4 || isempty(detectorFcn)
    detectorFcn = @computeMSLEdgeAndCore;
end

assert(isfile(labelFile), 'scoreMSLEdgeAgainstLabels:noLabels', ...
    'No label file found at %s; create one with labelMSLTopEdgeGUI.m first.', labelFile);
S = load(labelFile, 'labels');
labels = S.labels;
if ~isfield(labels, 'absentIntervals')
    for i = 1:numel(labels)
        labels(i).absentIntervals = datetime.empty(0,2); %#ok<AGROW>
    end
end

nTraced = sum(~cellfun(@isempty, {labels.pointsTime}));
nAbsent = sum([labels.isAbsent]);
nPartial = sum(cellfun(@(x) size(x,1)>0, {labels.absentIntervals}));
nUnlabeled = numel(labels) - nTraced - nAbsent;
fprintf(['scoreMSLEdgeAgainstLabels: %d traced (%d with a partial absent band), %d marked fully absent, ' ...
    '%d not yet labeled (of %d total).\n'], nTraced, nPartial, nAbsent, nUnlabeled, numel(labels));
assert(nTraced + nAbsent > 0, 'scoreMSLEdgeAgainstLabels:noLabeledWindows', ...
    'No windows have been labeled yet; label some with labelMSLTopEdgeGUI.m first.');

edge = detectorFcn(nascStitched, cfg);

%% Traced windows: depth error + miss rate (absent-interval time excluded)
perWinRows = {};
allErr = []; allMissed = 0; allPresentSteps = 0;
for i = 1:numel(labels)
    lbl = labels(i);
    if isempty(lbl.pointsTime) || lbl.isAbsent
        continue
    end
    inWin = edge.time >= lbl.winStart & edge.time <= lbl.winEnd & edge.deploymentID == lbl.deploymentID;
    tDet = edge.time(inWin);
    dDet = edge.edgeDepthTop(inWin);
    if numel(lbl.pointsTime) < 2
        truthAtDet = repmat(lbl.pointsDepth(1), size(tDet));
    else
        truthAtDet = interp1(lbl.pointsTime, lbl.pointsDepth, tDet, 'linear');
    end
    inLabelSpan = tDet >= min(lbl.pointsTime) & tDet <= max(lbl.pointsTime);
    inLabelSpan = inLabelSpan & ~inAnyInterval(tDet, lbl.absentIntervals);
    truthAtDet(~inLabelSpan) = nan;

    isMissed = inLabelSpan & isnan(dDet);
    isMatched = inLabelSpan & isfinite(dDet);
    errM = dDet(isMatched) - truthAtDet(isMatched);   % + = detector too deep

    perWinRows(end+1, :) = {lbl.deploymentID, lbl.winStart, lbl.winEnd, numel(lbl.pointsTime), ...
        sum(inLabelSpan), sum(isMatched), sum(isMissed), mean(abs(errM), 'omitnan'), ...
        mean(errM, 'omitnan'), sqrt(mean(errM.^2, 'omitnan'))}; %#ok<AGROW>

    allErr = [allErr; errM]; %#ok<AGROW>
    allMissed = allMissed + sum(isMissed);
    allPresentSteps = allPresentSteps + sum(inLabelSpan);
end
perWindow = cell2table(perWinRows, 'VariableNames', ...
    {'deploymentID','winStart','winEnd','nLabelPoints','nDetectorSteps','nMatched','nMissed','meanAbsErrorM','biasM','rmsErrorM'});

%% Absent intervals: false-positive rate
absentRows = {};
for i = 1:numel(labels)
    lbl = labels(i);
    ivals = getAbsentIntervalsForLabel(lbl);
    if isempty(ivals)
        continue
    end
    inDeploy = edge.deploymentID == lbl.deploymentID;
    for k = 1:size(ivals, 1)
        inIval = inDeploy & edge.time >= ivals(k,1) & edge.time <= ivals(k,2);
        dDet = edge.edgeDepthTop(inIval);
        nFP = sum(isfinite(dDet)); nS = numel(dDet);
        absentRows(end+1, :) = {lbl.deploymentID, lbl.winStart, lbl.winEnd, ivals(k,1), ivals(k,2), ...
            lbl.isAbsent, nFP, nS}; %#ok<AGROW>
    end
end
varNames = {'deploymentID','winStart','winEnd','intervalStart','intervalEnd','isWholeWindow','nFalsePositive','nSteps'};
if isempty(absentRows)
    % a label set with no "mark absent" windows or bands (every window fully
    % traced) gives a valid zero-row table, built directly

    absentIntervals = table('Size', [0, numel(varNames)], ...
        'VariableTypes', {'string','datetime','datetime','datetime','datetime','logical','double','double'}, ...
        'VariableNames', varNames);
else
    absentIntervals = cell2table(absentRows, 'VariableNames', varNames);
end

allFP = sum(absentIntervals.nFalsePositive);
allAbsentSteps = sum(absentIntervals.nSteps);

score.perWindow = perWindow;
score.absentIntervals = absentIntervals;
score.overallMeanAbsErrorM = mean(abs(allErr), 'omitnan');
score.overallBiasM = mean(allErr, 'omitnan');
score.overallRmsErrorM = sqrt(mean(allErr.^2, 'omitnan'));
score.overallMissRate = allMissed / max(allPresentSteps, 1);
score.overallFalsePositiveRate = allFP / max(allAbsentSteps, 1);
score.nTracedWindows = height(perWindow);
score.nAbsentIntervals = height(absentIntervals);

fprintf(['\n=== MSL top-edge detector score ===\n' ...
    'Mean |error| (m)       : %6.2f\n' ...
    'Bias (m, + = too deep) : %+6.2f\n' ...
    'RMS error (m)          : %6.2f\n' ...
    'Miss rate (%%)          : %6.1f\n' ...
    'False-positive rate (%%): %6.1f\n' ...
    '(%d traced windows, %d absent intervals)\n'], ...
    score.overallMeanAbsErrorM, score.overallBiasM, score.overallRmsErrorM, ...
    100*score.overallMissRate, 100*score.overallFalsePositiveRate, ...
    score.nTracedWindows, score.nAbsentIntervals);

if nUnlabeled > 0
    fprintf('\n(%d window(s) not yet labeled.)\n', nUnlabeled);
end
end

function ivals = getAbsentIntervalsForLabel(lbl)
if lbl.isAbsent
    ivals = [lbl.winStart, lbl.winEnd];
elseif ~isempty(lbl.absentIntervals)
    ivals = lbl.absentIntervals;
else
    ivals = datetime.empty(0,2);
end
end

function tf = inAnyInterval(t, intervals)
tf = false(size(t));
for k = 1:size(intervals, 1)
    tf = tf | (t >= intervals(k,1) & t <= intervals(k,2));
end
end
