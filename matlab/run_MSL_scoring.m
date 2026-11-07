%% run_MSL_scoring.m
% Scores the seed-and-walk detector (computeMSLEdgeAndCore.m) against the hand
% labels made with run_MSL_labeling.m. Re-run after changing the settings in
% mslDetectorConfig.m to test a parameter change.
%
% Reported: mean absolute error, bias (positive = detector too deep) and RMS
% error of the edge depth, the miss rate inside traced windows, and the
% false-positive rate inside spans marked absent.
%
% Inputs  (CCE1_DATA_DIR, default <repo>/data)
%   WBAT_export_5m/CCE1_WBAT_01_*.csv, CCE1_WBAT_hourly.mat
%   MSL/MSL_TopEdge_HandLabels.mat in CCE1_RESULTS_DIR (CCE1_LABEL_FILE
%   overrides this path)
% Output  (CCE1_RESULTS_DIR, default <repo>/results)
%   MSL/MSL_score_perWindow.csv, MSL/MSL_score_absentIntervals.csv,
%   MSL/MSL_score.mat

repoRoot = fileparts(fileparts(mfilename('fullpath')));
dataDir = getenv('CCE1_DATA_DIR');
if isempty(dataDir), dataDir = fullfile(repoRoot, 'data'); end
resultsDir = getenv('CCE1_RESULTS_DIR');
if isempty(resultsDir), resultsDir = fullfile(repoRoot, 'results'); end
outDir = fullfile(resultsDir, 'MSL');
if ~isfolder(outDir), mkdir(outDir); end
addpath(fileparts(mfilename('fullpath')));

nascStitched = loadMSLNightGrid(dataDir);

cfg = mslDetectorConfig();
cfg.edgeLabelFile = getenv('CCE1_LABEL_FILE');
if isempty(cfg.edgeLabelFile)
    cfg.edgeLabelFile = fullfile(outDir, 'MSL_TopEdge_HandLabels.mat');
end

score = scoreMSLEdgeAgainstLabels(nascStitched, cfg);

% window times are in the stitched pseudo-time; add real calendar dates for display
score.perWindow.realWinStart = toRealTime(score.perWindow.winStart, nascStitched);
score.perWindow.realWinEnd = toRealTime(score.perWindow.winEnd, nascStitched);
score.absentIntervals.realIntervalStart = toRealTime(score.absentIntervals.intervalStart, nascStitched);
score.absentIntervals.realIntervalEnd = toRealTime(score.absentIntervals.intervalEnd, nascStitched);

disp(score.perWindow(:, {'realWinStart','realWinEnd','nLabelPoints','nMatched','nMissed','meanAbsErrorM','biasM','rmsErrorM'}))
disp(score.absentIntervals(:, {'realIntervalStart','realIntervalEnd','nFalsePositive','nSteps'}))

save(fullfile(outDir, 'MSL_score.mat'), 'score', 'cfg');
writetable(score.perWindow, fullfile(outDir, 'MSL_score_perWindow.csv'));
writetable(score.absentIntervals, fullfile(outDir, 'MSL_score_absentIntervals.csv'));
fprintf('Saved score and per-window / absent-interval tables to %s\n', outDir);

function realT = toRealTime(stitchedT, nascStitched)
realT = NaT(size(stitchedT));
for k = 1:numel(stitchedT)
    [~, idx] = min(abs(nascStitched.time - stitchedT(k)));
    realT(k) = nascStitched.realTime(idx);
end
end
