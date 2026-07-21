function twix_hdr = read_twix_hdr_wrapper(inFile)
% wrapper to call read_twix_hdr which requires a file id already...
% gallichand@cardiff.ac.uk


fid = fopen(inFile,'r','l','US-ASCII');

% start of actual measurement data (sans header)
fseek(fid,0,'bof');

firstInt  = fread(fid,1,'uint32');
secondInt = fread(fid,1,'uint32');

% lazy software version check (VB or VD?)
if and(firstInt < 10000, secondInt <= 64)
    version = 'vd';
    disp('Software version: VD (!?)');

    % number of different scans in file stored in 2nd in
    NScans = secondInt;
    measID = fread(fid,1,'uint32');
    fileID = fread(fid,1,'uint32');
    % measOffset: points to beginning of header, usually at 10240 bytes
    measOffset = fread(fid,1,'uint64');
    measLength = fread(fid,1,'uint64');
    fseek(fid,measOffset,'bof');
    hdrLength  = fread(fid,1,'uint32');

else
    % in VB versions, the first 4 bytes indicate the beginning of the
    % raw data part of the file
    version  = 'vb';
    disp('Software version: VB (!?)');
    measOffset = 0;
    hdrLength  = firstInt;
    NScans     = 1; % VB does not support multiple scans in one file
end
datStart = measOffset + hdrLength;

%SRY read data correction factors
% do this for all VB datasets, so that the factors are available later
% in the image_obj if the user chooses to set the correction flag
if (strcmp(version, 'vb')) % not implemented/tested for vd, yet
    frewind(fid);
    while ( (ftell(fid) < datStart) && ~exist('rawfactors', 'var'))
        line = fgetl(fid);
        %find the section of the protocol
        %note: the factors are also available in <ParamArray."CoilSelects">
        %along with element name and FFT scale, but the parsing is
        %significantly more difficult
        if (~isempty(strfind(line, '<ParamArray."axRawDataCorrectionFactor">')))
            while (ftell(fid) < datStart)
                line = fgetl(fid);
                %find the line with correction factors
                %the factors are on the first line that begins with this
                %pattern
                if (~isempty(strfind(line, '{ {  { ')))
                    line = strrep(line, '}  { } ', '0.0');
                    line = strrep(line, '{', '');
                    line = strrep(line, '}', '');
                    rawfactors = textscan(line, '%f');
                    rawfactors = rawfactors{1}; %textscan returns a 1x1 cell array
                    % this does not work in this location because the MDHs
                    % have not been parsed yet
                    %                    if (length(rawfactors) ~= 2*max(image_obj.NCha))
                    %                       error('Number of raw factors (%f) does not equal channel count (%d)', length(rawfactors)/2, image_obj.NCha);
                    %                    end;
                    if (mod(length(rawfactors),2) ~= 0)
                        error('Error reading rawfactors');
                    end;
                    %note the transpose, this makes the vector
                    %multiplication during the read easier
                    arg.rawDataCorrectionFactors = rawfactors(1:2:end).' + 1i*rawfactors(2:2:end).';
                    break;
                end
            end
        end
    end
    disp('Read raw data correction factors');
end

% data will be read in two steps (two while loops):
%   1) reading all MDHs to find maximum line no., partition no.,... for
%      ima, ref,... scan
%   2) reading the data
tic;
cPos            = measOffset;
twix_hdr        = cell(1,NScans);

for s=1:NScans
    fseek(fid,cPos,'bof');
    hdr_len = fread(fid, 1,'uint32');

    [twix_hdr{s},rstraj] = read_twix_hdr(fid);
    
end
if NScans == 1
    twix_hdr = twix_hdr{1};
end
fclose(fid);