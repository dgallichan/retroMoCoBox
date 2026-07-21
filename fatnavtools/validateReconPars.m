function [isValid, errMsg] = validateReconPars(theseReconPars)
% function validateReconPars(theseReconPars)
%
% Check that the current parameters have all the required fields

defaultReconPars = getDefaultReconPars();

isValid = 1;

onlyInDefault = setdiff(fieldnames(defaultReconPars),fieldnames(theseReconPars));

if ~isempty(onlyInDefault)
    isValid = 0;
    onlyInDefault = strjoin(onlyInDefault, newline + "  - ");
    errMsg1 = sprintf('Missing parameters:\n  - %s \n',onlyInDefault);
else
    errMsg1 = "";
end

onlyInThese = setdiff(fieldnames(theseReconPars),fieldnames(defaultReconPars));
if ~isempty(onlyInThese)
    isValid = 0;
    onlyInThese = strjoin(onlyInThese, newline + "  - ");
    errMsg2 = sprintf('Unknown parameters: \n  - %s \n',onlyInThese);
else
    errMsg2 = "";
end

if ~isValid
    errMsg = "" + newline + newline + "**** Error: reconPars struct does not appear valid ****" + newline + newline;
    errMsg = errMsg + errMsg1 + errMsg2;
    errMsg = errMsg + newline + "Did you generate defaults with reconPars = getDefaultReconPars() ?";
else
    errMsg = "";
end
