function MIDstr = getMIDstr(thisFile)

iUnderscore = strfind(thisFile,'_');
iMID = strfind(thisFile,'_MID');
iUnderscore(iUnderscore<=iMID) = [];
MIDstr = thisFile(iMID+1:iUnderscore(1)-1);