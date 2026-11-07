function [nascStitched, night] = loadMSLNightGrid(dataDir)
%LOADMSLNIGHTGRID Night-time NASC grid for the MSL tracing, labeling and scoring.
% Reads the monthly 5 m exports in dataDir/WBAT_export_5m and the night-time
% windows (variable night_all) in dataDir/CCE1_WBAT_hourly.mat, then keeps the
% night-time rows trimmed by 1 h at both ends and stitched into one series.
%
% OUTPUTS
%   nascStitched  struct from restrictToNightTrimmed.m (time, realTime,
%                 nightIndex, depth, NASC)
%   night         [nNights x 2] datenum, sunset and following sunrise

nasc = buildNASC5mGrid(fullfile(dataDir, 'WBAT_export_5m'));
S = load(fullfile(dataDir, 'CCE1_WBAT_hourly.mat'), 'night_all');
night = S.night_all;
nascStitched = restrictToNightTrimmed(nasc, night, 1);
end
