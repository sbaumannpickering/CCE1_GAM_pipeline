%% run_MSL_seedwalk.m
% Traces the top edge, core depth and integrated backscatter of the
% mesopelagic scattering layer (MSL) for every night of the CCE1 record and
% writes MSL_night_seedWalk.csv, the input of R/08_watermass_regression.R.
%
% Steps
%   1. Load the monthly 5 m WBAT exports into one hourly NASC x depth grid
%      (buildNASC5mGrid.m).
%   2. Keep the night-time rows, trimmed by 1 h at both ends of every night
%      and stitched into one continuous series (restrictToNightTrimmed.m).
%   3. Seed-and-walk detector for the layer top edge and core
%      (computeMSLEdgeAndCore.m).
%   4. Per-night median of the hours with a detection.
%   5. Gap filling of short runs (<= 2 nights, not touching either end of the
%      record; interpolateShortGaps.m).
%   6. Five-night running median of the gap-filled series.
%
% Inputs  (CCE1_DATA_DIR, default <repo>/data)
%   WBAT_export_5m/CCE1_WBAT_01_*.csv
%   CCE1_WBAT_hourly.mat  (variable night_all: night-time windows)
% Output  (CCE1_RESULTS_DIR, default <repo>/results)
%   MSL/MSL_night_seedWalk.csv
%
% Needs MATLAB (written for R2025b) and the Statistics and Machine Learning
% Toolbox (prctile, mad).

repoRoot = fileparts(fileparts(mfilename('fullpath')));
dataDir = getenv('CCE1_DATA_DIR');
if isempty(dataDir), dataDir = fullfile(repoRoot, 'data'); end
resultsDir = getenv('CCE1_RESULTS_DIR');
if isempty(resultsDir), resultsDir = fullfile(repoRoot, 'results'); end
outDir = fullfile(resultsDir, 'MSL');
if ~isfolder(outDir), mkdir(outDir); end
addpath(fileparts(mfilename('fullpath')));

%% 1-2. Night-time grid
[nascStitched, night] = loadMSLNightGrid(dataDir);

%% 3. Detector
edge = computeMSLEdgeAndCore(nascStitched, mslDetectorConfig());

%% 4. Per-night aggregation (median of each night's hours with a detection)
nNights = size(night, 1);
nightStart = datetime(night(:,1), 'ConvertFrom', 'datenum');
nightEnd = datetime(night(:,2), 'ConvertFrom', 'datenum');
medianEdgeM = nan(nNights,1); medianCoreDepthM = nan(nNights,1);
medianCoreAmplitude = nan(nNights,1); medianIntegratedNASC = nan(nNights,1);
nHoursPresent = zeros(nNights,1); nHoursTotal = zeros(nNights,1);

for n = 1:nNights
    inNight = edge.nightIndex == n;
    nHoursTotal(n) = sum(inNight);
    edgeVals = edge.edgeDepthTop(inNight);
    coreDepthVals = edge.coreDepth(inNight);
    coreAmpVals = edge.coreAmplitude(inNight);
    nascVals = edge.integratedNASC(inNight);
    present = isfinite(edgeVals);
    nHoursPresent(n) = sum(present);
    if any(present)
        medianEdgeM(n) = median(edgeVals(present));
        medianCoreDepthM(n) = median(coreDepthVals(present));
        medianCoreAmplitude(n) = median(coreAmpVals(present));
        medianIntegratedNASC(n) = median(nascVals(present));
    end
end
pctPresent = 100 * nHoursPresent ./ max(nHoursTotal, 1);

%% 5. Fill short gaps (<= 2 nights). Real sunset times jitter by minutes from
% night to night, which would fail the regular-spacing check of
% interpolateShortGaps.m, so a synthetic axis of one day per night is used
% (the night index is what is regular).
syntheticAxis = datetime(2000,1,1) + days(0:nNights-1)';
maxGapNights = days(2);

[edgeFilled, edgeWasInterp] = interpolateShortGaps(syntheticAxis, medianEdgeM, maxGapNights);
[coreDepthFilled, ~] = interpolateShortGaps(syntheticAxis, medianCoreDepthM, maxGapNights);
[coreAmpFilled, ~] = interpolateShortGaps(syntheticAxis, medianCoreAmplitude, maxGapNights);
[integNASCFilled, ~] = interpolateShortGaps(syntheticAxis, medianIntegratedNASC, maxGapNights);

%% 6. Five-night persistence filter on the gap-filled series (nights still
% missing after gap filling stay missing)
persistWindowNights = 5;
wasPresentAfterFill = isfinite(edgeFilled);
filteredEdgeM = movmedian(edgeFilled, persistWindowNights, 'omitnan'); filteredEdgeM(~wasPresentAfterFill) = nan;
filteredCoreDepthM = movmedian(coreDepthFilled, persistWindowNights, 'omitnan'); filteredCoreDepthM(~wasPresentAfterFill) = nan;
filteredCoreAmplitude = movmedian(coreAmpFilled, persistWindowNights, 'omitnan'); filteredCoreAmplitude(~wasPresentAfterFill) = nan;
filteredIntegratedNASC = movmedian(integNASCFilled, persistWindowNights, 'omitnan'); filteredIntegratedNASC(~wasPresentAfterFill) = nan;

nightTable = table(nightStart, nightEnd, medianEdgeM, edgeWasInterp, edgeFilled, filteredEdgeM, ...
    medianCoreDepthM, coreDepthFilled, filteredCoreDepthM, ...
    medianCoreAmplitude, coreAmpFilled, filteredCoreAmplitude, ...
    medianIntegratedNASC, integNASCFilled, filteredIntegratedNASC, ...
    nHoursPresent, nHoursTotal, pctPresent);
outCsv = fullfile(outDir, 'MSL_night_seedWalk.csv');
writetable(nightTable, outCsv);
fprintf('Saved %s\n', outCsv);
fprintf('%d/%d nights with a raw detection (%.1f%%); %d/%d after short-gap fill (%.1f%%).\n', ...
    sum(isfinite(medianEdgeM)), nNights, 100*sum(isfinite(medianEdgeM))/nNights, ...
    sum(wasPresentAfterFill), nNights, 100*sum(wasPresentAfterFill)/nNights);
