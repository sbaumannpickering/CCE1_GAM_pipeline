%% run_MSL_labeling.m
% Opens the hand-labeling tool (labelMSLTopEdgeGUI.m) on the night-time grid, to
% trace where you judge the layer top edge to be. The labels are the ground truth
% for run_MSL_scoring.m. Run it in an interactive MATLAB session (the tool needs
% a window; matlab -batch does not work).
%
% 16 windows of 50 stitched hours (about 5 nights each) are spread evenly across
% the record. Progress is saved on every navigation and when the window is
% closed; running the script again resumes at the first unlabeled window.
% Changing the window settings later keeps labels whose window start time still
% matches and backs up the old file.
%
% Inputs  (CCE1_DATA_DIR, default <repo>/data)
%   WBAT_export_5m/CCE1_WBAT_01_*.csv, CCE1_WBAT_hourly.mat
% Output  (CCE1_RESULTS_DIR, default <repo>/results)
%   MSL/MSL_TopEdge_HandLabels.mat
%   (CCE1_LABEL_FILE overrides this path)

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
cfg.edgeLabelWindowHours = 50;
cfg.edgeLabelWindowsPerDeployment = 16;
cfg.edgeLabelFile = getenv('CCE1_LABEL_FILE');
if isempty(cfg.edgeLabelFile)
    cfg.edgeLabelFile = fullfile(outDir, 'MSL_TopEdge_HandLabels.mat');
end

labelMSLTopEdgeGUI(nascStitched, cfg);
