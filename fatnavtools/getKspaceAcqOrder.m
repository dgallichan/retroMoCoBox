function [idxMap_acrossEchoTrains, idxMap_withinEchoTrain, lineCounterMap, timeMap, rawCorrectMap] = getKspaceAcqOrder(twix_obj)
%
% This function currently appears to work for various SPACE acquisitions
% with sagittal acquistion and readout in AP dirn

idxMap_acrossEchoTrains = zeros(twix_obj.image.NLin,twix_obj.image.NPar);
idxMap_withinEchoTrain = zeros(twix_obj.image.NLin,twix_obj.image.NPar);
timeMap = zeros(twix_obj.image.NLin,twix_obj.image.NPar);
rawCorrectMap = zeros(twix_obj.image.NLin,twix_obj.image.NPar);
lineCounterMap = zeros(twix_obj.image.NLin,twix_obj.image.NPar);


for iAcq = 1:twix_obj.imageWithRefscan.NAcq
    idxMap_withinEchoTrain(twix_obj.imageWithRefscan.Lin(iAcq),twix_obj.imageWithRefscan.Par(iAcq)) = ...
        twix_obj.imageWithRefscan.Seg(iAcq);

    timeMap(twix_obj.imageWithRefscan.Lin(iAcq),twix_obj.imageWithRefscan.Par(iAcq)) = ...
        twix_obj.imageWithRefscan.timestamp(iAcq);

    rawCorrectMap(twix_obj.imageWithRefscan.Lin(iAcq),twix_obj.imageWithRefscan.Par(iAcq)) = ...
        twix_obj.imageWithRefscan.IsRawDataCorrect(iAcq);
    
    lineCounterMap(twix_obj.imageWithRefscan.Lin(iAcq),twix_obj.imageWithRefscan.Par(iAcq)) = ...
        lineCounterMap(twix_obj.imageWithRefscan.Lin(iAcq),twix_obj.imageWithRefscan.Par(iAcq)) +1;
end


% This method here looks questionable, as NAcq doesn't divide by NSeg
% evenly!
nPerSeg = ceil(twix_obj.imageWithRefscan.NAcq / twix_obj.imageWithRefscan.NSeg);
for iSeg = 1:nPerSeg
    thisRange = [1:twix_obj.imageWithRefscan.NSeg]+(iSeg-1)*twix_obj.imageWithRefscan.NSeg;
    thisRange(thisRange>twix_obj.imageWithRefscan.NAcq) = [];
    for iR = 1:length(thisRange)
        idxMap_acrossEchoTrains(twix_obj.imageWithRefscan.Lin(thisRange(iR)),...
            twix_obj.imageWithRefscan.Par(thisRange(iR))) = iSeg;
    end
end


% % this method seems more 'correct' (but not necessarily robust!) - but in
% % the test dataset gave the same answer as above - it looks like it is
% % only the final train in the entire acquisition which is 'cut short'
% newIdxMap = zeros(twix_obj.image.NLin,twix_obj.image.NPar);
% iTrain = 1;
% for iAcq = 1:twix_obj.imageWithRefscan.NAcq
%     if iAcq > 1
%         timeSinceLastAcq = twix_obj.imageWithRefscan.timestamp(iAcq) - twix_obj.imageWithRefscan.timestamp(iAcq-1);
%     else 
%         timeSinceLastAcq = 0;
%     end
%     
%     if timeSinceLastAcq > 100 % arbitrary here, looks from testing like I get 2 or 3 for in-train and 1699/1700 for across
%         iTrain = iTrain + 1;
%     end
%     
%     newIdxMap(twix_obj.imageWithRefscan.Lin(iAcq),twix_obj.imageWithRefscan.Par(iAcq)) = iTrain;
%     
% end
% newIdxMap = newIdxMap.';


%%% can transpose here if preferred:
% idxMap_acrossEchoTrains = idxMap_acrossEchoTrains.';
% idxMap_withinEchoTrain = idxMap_withinEchoTrain.';
% lineCounterMap = lineCounterMap.';
% timeMap = timeMap.';
% rawCorrectMap = rawCorrectMap.';

timeMap(timeMap>0) = timeMap(timeMap>0) - min(timeMap(timeMap>0)) + 1; % reset to more sensible numbers

if any(lineCounterMap(:)>1)
    warning('Seems some lines were acquired more than once!')
end
