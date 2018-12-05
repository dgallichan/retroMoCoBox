function MIDstr = getMIDstr(thisFile)
% function MIDstr = getMIDstr(thisFile)
%
% From a filename of the standard Siemens TWIX form, e.g.
% meas_MID00100_FID03307_MPRAGE_motion_FatNavs_1mm_iso.dat
% this function should return the number after the letters MID. This is
% useful to keep track of separate datasets in the same folder, and is used
% to name the output folder and motion parameters.
% 
% If the file doesn't contain '_MID' then MIDstr will return as simply 'XX'

iUnderscore = strfind(thisFile,'_');
iMID = strfind(thisFile,'_MID');
if isempty(iMID)
    MIDstr = 'XX';
else
    iUnderscore(iUnderscore<=iMID) = [];
    MIDstr = thisFile(iMID+1:iUnderscore(1)-1);
end