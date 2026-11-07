function [filled, wasInterpolated] = interpolateShortGaps(time, x, maxGapDuration)
%INTERPOLATESHORTGAPS Linearly interpolate short NaN runs in a regularly
% spaced series. NaN runs no longer than maxGapDuration are filled; longer
% runs, and any run touching either end of the record (which cannot be
% bracketed), are left as NaN. time must be regularly spaced (asserted; the
% median step is the sampling interval). Interpolation is by sample index,
% which is linear in time for a regular series.
%
% INPUTS
%   time           [n x 1] datetime, regularly spaced
%   x              [n x 1] double, may contain NaN runs
%   maxGapDuration duration scalar (e.g. days(2))
%
% OUTPUTS
%   filled          [n x 1] double, x with short NaN runs filled
%   wasInterpolated [n x 1] logical, true where filled(i) was interpolated
%                   (false for original finite values and for NaN runs left
%                   unfilled)
time = time(:);
x = x(:);
n = numel(x);
assert(numel(time) == n, 'interpolateShortGaps:sizeMismatch', 'time and x must be the same length.');
assert(n >= 2, 'interpolateShortGaps:tooFewPoints', 'Need at least 2 samples.');

dt = median(diff(time));
stepOK = abs(diff(time) - dt) < dt * 0.01;
assert(all(stepOK), 'interpolateShortGaps:notRegular', ...
    'time must be regularly spaced (median step %s); found irregular step(s) at %d/%d transitions.', ...
    string(dt), sum(~stepOK), numel(stepOK));

isNaNx = isnan(x);
filled = x;
wasInterpolated = false(n, 1);

if ~any(isNaNx)
    fprintf('interpolateShortGaps: no NaN values present -- nothing to do.\n');
    return
end

padded = [false; isNaNx; false];
d = diff(padded);
runStarts = find(d == 1);
runEnds   = find(d == -1) - 1;
runDurationSamples = runEnds - runStarts + 1;
runDuration = runDurationSamples * dt;

nFilled = 0; nLeftShort = 0; nLeftLong = 0;
samplesFilled = 0; samplesLeftShort = 0; samplesLeftLong = 0;
for r = 1:numel(runStarts)
    idx = runStarts(r):runEnds(r);
    touchesEdge = runStarts(r) == 1 || runEnds(r) == n;
    if touchesEdge
        % Can't bracket a run with no finite value on one side, regardless
        % of duration -- left as NaN.
        nLeftShort = nLeftShort + 1; samplesLeftShort = samplesLeftShort + numel(idx);
        continue
    end
    if runDuration(r) <= maxGapDuration
        i0 = runStarts(r) - 1; i1 = runEnds(r) + 1;
        x0 = x(i0); x1 = x(i1);
        frac = (idx - i0) / (i1 - i0);
        filled(idx) = x0 + (x1 - x0) .* frac;
        wasInterpolated(idx) = true;
        nFilled = nFilled + 1; samplesFilled = samplesFilled + numel(idx);
    else
        nLeftLong = nLeftLong + 1; samplesLeftLong = samplesLeftLong + numel(idx);
    end
end

fprintf(['interpolateShortGaps: %d NaN run(s) found (dt=%s). %d run(s), %d sample(s) interpolated (<= %s, ' ...
    'not touching a record edge). %d run(s), %d sample(s) left as NaN (> %s). %d run(s), %d sample(s) left as ' ...
    'NaN (touch a record edge, unbracketable regardless of duration).\n'], numel(runStarts), string(dt), ...
    nFilled, samplesFilled, string(maxGapDuration), nLeftLong, samplesLeftLong, string(maxGapDuration), ...
    nLeftShort, samplesLeftShort);
end
