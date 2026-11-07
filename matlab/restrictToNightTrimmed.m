function nascStitched = restrictToNightTrimmed(nasc, night, trimHours)
%RESTRICTTONIGHTTRIMMED Keep only night-time rows of a NASC grid, trimmed at
% both ends of every night, and stitch them into one continuous time axis.
%
% trimHours is removed from BOTH ends of every sunset-sunrise interval,
% because the first and last hour of night still show active migration, which
% can look like a false layer edge. The kept rows are stitched into one evenly
% spaced pseudo-time axis (nasc.time) with no daytime gap, so the detector's
% elapsed-time logic sees one ordinary hourly series and can carry a track
% across nights. The true calendar time of each kept row is returned in
% .realTime, and .nightIndex gives the night (1..nNights) each row belongs to.
%
% INPUTS
%   nasc       struct with time [nTime x 1 datetime], depth [1 x nDepth] and
%              NASC [nTime x nDepth] (optional deploymentID)
%   night      [nNights x 2] night-time windows: sunset in column 1 and the
%              following sunrise in column 2 (MATLAB datenum). The variable
%              `night_all` in CCE1_WBAT_hourly.mat.
%   trimHours  hours to trim off each end of every night (1 is used)
%
% OUTPUT  struct with time (stitched), realTime, nightIndex, depth, NASC and,
%         if present, deploymentID.

nNights = size(night, 1);

timeDatenum = datenum(nasc.time); %#ok<DATNM>
nightIndexAll = zeros(size(timeDatenum));
for n = 1:nNights
    lo = night(n,1) + trimHours/24;
    hi = night(n,2) - trimHours/24;
    inNight = timeDatenum >= lo & timeDatenum <= hi;
    nightIndexAll(inNight) = n;
end

keep = nightIndexAll > 0;
realTime = nasc.time(keep);
nightIndex = nightIndexAll(keep);
NASC = nasc.NASC(keep, :);

% realTime is already chronologically sorted (nasc.time was ascending and
% keep is a boolean mask on it), so a uniform per-kept-row increment
% stitches them into one continuous series with no information about the
% true gap sizes -- deliberately: that is what "without a break" means
% here, treating elapsed-sample-count as the detector's only notion of
% "time passed" rather than true wall-clock hours.
nKept = sum(keep);
stitchedTime = realTime(1) + hours(0:nKept-1)';

nascStitched = struct();
nascStitched.time = stitchedTime;
nascStitched.realTime = realTime;
nascStitched.nightIndex = nightIndex;
nascStitched.depth = nasc.depth;
nascStitched.NASC = NASC;
if isfield(nasc, 'deploymentID') && ~isempty(nasc.deploymentID)
    nascStitched.deploymentID = nasc.deploymentID(keep);
end

fprintf(['restrictToNightTrimmed: kept %d/%d steps (%.1f%%) after trimming %.1f h off each end of %d nights; ' ...
    'stitched into one continuous %d-step series with no day gaps.\n'], ...
    nKept, numel(timeDatenum), 100*nKept/numel(timeDatenum), trimHours, nNights, nKept);
end
