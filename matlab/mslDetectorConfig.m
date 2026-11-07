function cfg = mslDetectorConfig()
%MSLDETECTORCONFIG Detector settings for computeMSLEdgeAndCore.m, shared by
% the tracing, labeling and scoring scripts. Other detector settings (layer
% thickness for smoothing, continuity gate, timeout, re-acquisition) use the
% defaults documented in computeMSLEdgeAndCore.m.

cfg.edgeSearchDepthRangeM = [100 495];   % m, where the layer edge and core are searched
cfg.edgeRefDepthRangeM    = [50 100];    % m, quiet reference water above the layer
cfg.edgeZHigh = 4;                       % z-score the core must clear
cfg.edgeZLow = 2.5;                      % z-score a cell must clear to extend the walk
cfg.edgeMinAbsNASCPercentile = 50;       % absolute floor (percentile of positive NASC)
cfg.edgeAbsFloorLowFrac = 0.85;          % fraction of the floor a cell must clear to extend
end
