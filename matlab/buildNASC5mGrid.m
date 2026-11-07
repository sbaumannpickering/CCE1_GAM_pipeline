function nasc = buildNASC5mGrid(dataDir, excludeLayerMin, excludeLayerMax)
%BUILDNASC5MGRID Load the monthly 5 m WBAT exports (CCE1_WBAT_01_*.csv) into
% one native-resolution NASC struct: time ascending, depth ascending
% (shallow to deep), linear NASC.
%
% Each row of an export is one (Interval, Layer) cell. Interval is one hourly
% time step and Layer is one 5 m depth bin. Depth_mean is metres below the
% surface (transducer depth 504 m minus range from the transducer), so
% increasing Layer means decreasing depth. Date_M and Time_M are repeated
% across all layers of an interval.
%
% Excluded layers (instrument artifacts at both ends of the beam):
%   Layer < excludeLayerMin  the bin at the transducer face (Layer -100,
%                            ringdown, NASC ~4e9 m^2 nmi^-2)
%   Layer > excludeLayerMax  bins past the sea surface (Layer 1 and above,
%                            surface reflection). The median NASC profile is
%                            smooth through Layer 0 and jumps by an order of
%                            magnitude at Layer 1 and again at Layer 2.
%
% INPUTS
%   dataDir          folder containing CCE1_WBAT_01_*.csv
%   excludeLayerMin  lowest Layer index kept (default -99)
%   excludeLayerMax  highest Layer index kept (default 0)
%
% OUTPUT  struct nasc
%   time   [nTime x 1] datetime, ascending, one per Interval present
%          (irregular spacing around deployment turnarounds)
%   depth  [1 x nDepth] double, metres below the surface, ascending
%   NASC   [nTime x nDepth] double, linear NASC; NaN where a cell is missing

if nargin < 2 || isempty(excludeLayerMin), excludeLayerMin = -99; end
if nargin < 3 || isempty(excludeLayerMax), excludeLayerMax = 0; end

files = dir(fullfile(dataDir, 'CCE1_WBAT_01_*.csv'));
assert(~isempty(files), 'buildNASC5mGrid:noFiles', 'No CCE1_WBAT_01_*.csv files found in %s', dataDir);

allTime = {};
allLayer = {};
allNASC = {};
layerToDepth = containers.Map('KeyType', 'int32', 'ValueType', 'double');

for f = 1:numel(files)
    fp = fullfile(files(f).folder, files(f).name);
    T = readtable(fp, 'TextType', 'string');

    keep = T.Layer >= excludeLayerMin & T.Layer <= excludeLayerMax;
    T = T(keep, :);

    % depth of each layer, recorded once (identical across files)
    uL = unique(T.Layer);
    for i = 1:numel(uL)
        L = int32(uL(i));
        if ~isKey(layerToDepth, L)
            layerToDepth(L) = T.Depth_mean(find(T.Layer == uL(i), 1));
        end
    end

    % per-row datetime from Date_M (yyyymmdd, numeric) + Time_M (HH:MM:SS.ffff)
    dateOnly = datetime(string(T.Date_M), 'InputFormat', 'yyyyMMdd');
    if isduration(T.Time_M)
        dt = dateOnly + T.Time_M;
    else
        dt = dateOnly + duration(strtrim(string(T.Time_M)), 'InputFormat', 'hh:mm:ss.SSSS');
    end

    allTime{end+1} = dt; %#ok<AGROW>
    allLayer{end+1} = T.Layer; %#ok<AGROW>
    allNASC{end+1} = T.NASC; %#ok<AGROW>

    fprintf('Loaded %s (%d/%d), %d rows.\n', files(f).name, f, numel(files), height(T));
end

timeCat = vertcat(allTime{:});
layerCat = vertcat(allLayer{:});
nascCat = vertcat(allNASC{:});

[uTime, ~, timeIdx] = unique(timeCat);
layerKeys = sort(cell2mat(keys(layerToDepth)));
nDepth = numel(layerKeys);
layerIdxMap = containers.Map(num2cell(layerKeys), num2cell(1:nDepth));

NASC = nan(numel(uTime), nDepth);
layerColIdx = cell2mat(values(layerIdxMap, num2cell(int32(layerCat))));
NASC(sub2ind(size(NASC), timeIdx, layerColIdx(:))) = nascCat;

% ascending Layer is descending depth: reorder the columns shallow to deep
depth = cell2mat(values(layerToDepth, num2cell(layerKeys)));
[depth, order] = sort(depth);
NASC = NASC(:, order);

nasc = struct();
nasc.time = uTime;
nasc.depth = depth(:)';
nasc.NASC = NASC;

fprintf('Final grid: %d time steps x %d depth bins (%.1f-%.1f m), %s to %s.\n', ...
    numel(nasc.time), nDepth, min(nasc.depth), max(nasc.depth), ...
    datestr(nasc.time(1)), datestr(nasc.time(end))); %#ok<DATST>
end
