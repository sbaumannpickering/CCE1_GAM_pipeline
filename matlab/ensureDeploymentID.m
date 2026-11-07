function nasc = ensureDeploymentID(nasc)
%ENSUREDEPLOYMENTID Give a NASC grid struct a deploymentID field if it has
% none -- every row gets "ALL", i.e. the whole record is treated as one
% continuous segment. The detector, label-window builders, GUI and scorer
% all key on deploymentID (tracking state is reset at every change in it),
% so a single-deployment dataset needs it filled in, not left out.
if ~isfield(nasc, 'deploymentID') || isempty(nasc.deploymentID)
    nasc.deploymentID = repmat("ALL", numel(nasc.time), 1);
else
    nasc.deploymentID = string(nasc.deploymentID(:));
end
end
