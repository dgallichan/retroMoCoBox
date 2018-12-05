function [timingReport, tempNameRoots] = performHostGRAPPArecon_toDisk(twix_obj,tempDir,counterStruct, nSliceNeighbours, combinePars)
%
% function [timingReport, tempNameRoots] = performHostGRAPPArecon_toDisk(twix_obj,tempDir,counterStruct, nSliceNeighbours, combinePars)
%
% Applies GRAPPA to the host data, which here is assumed to have acceleration in
% the first phase-encoding direction only (as is the current limitation of
% the MP2RAGE Siemens code).
% 
% %% Will only currently work with 1D GRAPPA...!
% %%%%%%%%%%%%%%%%%%%
%
% A 1D FFT is performed in the readout direction and a separate temporary file 
% is created for each virtual 'slice' in that direction - massively reducing RAM requirements
% (but significantly slower!).
%
% The _toDisk.m variation is supposed to be parfor capable, but still
% writes files to disk
%
% -------
% daniel.gallichan@epfl.ch, March 2016
%
% Update Jan 2018, gallichand@cardiff.ac.uk 
% -- clean up for actual integration into retroMoCoBox pipeline
%
% -- Feb 2018 - Adapt to handle 'averages' and 'repetitions'
%    - not that it isn't straightforward to keep track of filenames of the
%    temporary files to make sure nothing bad happens if you run this
%    simultaneously on multiple averages or repetitions - probably best if
%    they each have a separate 'tempDir'

if nargin < 5
    combinePars = [];
end

if nargin < 4
    nSliceNeighbours = 2; % no. of neighbouring 'virtual' slices to include in ACS lines to try to improve conditioning and regularize in the readout direction
    % quick tests suggest that 2 is a good number to use... :)
end

if nargin < 3
    counterStruct.iRep = 1;
    counterStruct.iAve = 1;
else
    if ~isfield(counterStruct,'iRep')
        counterStruct.iRep = 1;
    end
    if ~isfield(counterStruct,'iAve')
        counterStruct.iAve = 1;
    end
end



%%

nread = twix_obj.image.dataSize(1);
npe1 = twix_obj.image.dataSize(3);
npe2 = twix_obj.image.dataSize(4);
nc = twix_obj.image.dataSize(2);
nSet = twix_obj.image.dataSize(10);
nEco = twix_obj.image.NEco;

nAve = twix_obj.image.NAve;
nRep = twix_obj.image.NRep;

if counterStruct.iAve > nAve
    disp(['Error, tried to recon average ' num2str(counterStruct.iAve) ', but only found ' num2str(nAve) ' averages in data']);
    return
end
if counterStruct.iRep > nRep
    disp(['Error, tried to recon repetition ' num2str(counterStruct.iRep) ', but only found ' num2str(nRep) ' repetitions in data']);
    return
end


nACSmeas1 = twix_obj.refscan.dataSize(3);
nACSmeas2 = twix_obj.refscan.dataSize(4);

nACSmax = 40;
nACS1 = min(nACSmeas1,nACSmax);
nACS2 = min(nACSmeas2,nACSmax);

centerPar1 = twix_obj.image.centerPar(1);
centerLin1 = twix_obj.image.centerLin(1);

fileName = twix_obj.image.filename;
MIDstr = getMIDstr(fileName);


%%

if ~exist(tempDir,'dir')
    mkdir(tempDir)
end


%% make a separate temporary file per readout 'slice'

tempFiles = dir([tempDir '/temp_slicedata_*.bin']);

if length(tempFiles)~=nread % Check if we can assume they already exist...
    
    % open a file for each slice (max allowed number of simultaneous files on
    % my linux system seems to be 1024 ( ulimit -a) )
    fileid = zeros(nread,1);
    for iS = 1:nread
        fileid(iS) = fopen([tempDir '/temp_slicedata_' num2str(iS,'%.4d') '.bin'],'w');
    end
    
    % now go through the raw data and write out these files
    %%% this currently assumes that the whole data file divided by the size in
    %%% the second phase-encoding direction will fit
    %%% in RAM on your computer...
    disp('..............')
    disp('... Taking the raw data file, doing a 1D-FFT and then resaving a new file per ''slice''')
    tic
    for iPE2 = 1:npe2
        
        iPE1 = 1:npe1;
        iSet = 1:nSet;
        iEco = 1:nEco;
        
        thisdata = twix_obj.image(:,:,iPE1,iPE2,1,counterStruct.iAve,1,iEco,counterStruct.iRep,iSet);
        thisdata = ifft1s(thisdata,1);
        
        for iS = 1:nread
            fwrite(fileid(iS),real(thisdata(iS,:)),'single');
            fwrite(fileid(iS),imag(thisdata(iS,:)),'single');
        end
        if iPE2 > 1
            fprintf([repmat('\b',1,3)]);
        end
        fprintf('%s%%',num2str(round(100*iPE2/npe2),'%.2d'));
    end
    fprintf('\n');
    disp('Done')
    disp('..............')
    timingReport.FFTperSlice = toc; % 328s for 9 Gb file...
    
    % now close all those files
    for iS = 1:nread
        fclose(fileid(iS));
    end
    
end

%%% It would seem more elegant to also use a matfile variable for this stage, but it just seems to be *really* slow compared to the method above...

% mOut.data1DFFT = complex(zeros(nread, nc, npe1, npe2, nSet,'single'));
%
% tic
% for iPE2 = 1:npe2;
%     iPE1 = 1:npe1;
%     iSet = 1:nSet;
%
%     thisdata = twix_obj.image(:,:,iPE1,iPE2,1,1,1,1,1,iSet);
%     thisdata = squeeze(ifft1s(thisdata,1));
%
%     if nSet > 1
%         mOut.data1DFFT(:,:,iPE1,iPE2,iSet) = reshape(thisdata,[nread nc npe1 1 nSet]);
%     else
%         mOut.data1DFFT(:,:,iPE1,iPE2) =  thisdata;
%     end
%         fprintf('.')
% end
% fprintf('\n')
% disp('Done')
% disp('..............')
% timingReport.FFTperSlice = toc;



%% Now do GRAPPA on each of those slices


ACSdata = twix_obj.refscan(:,:,:,:,:,counterStruct.iAve,:,:,counterStruct.iRep,:);
ACSdata = ifft1s(ACSdata);

%%% Define GRAPPA kernel
gx = 4; gy = 3; % size of GRAPPA kernel
Ax = twix_obj.hdr.MeasYaps.sPat.lAccelFactPE;
Ay = 1;
grapKernel = zeros(gx + (gx-1)*(Ax-1), gy + (gy-1) * (Ay-1));

% source points are marked with 1
grapKernel(1:Ax:end, 1:Ay:end) = 1;

% target points are marked with 0.5
if gx == 1,    startx = 2; else  startx = 2 + (floor(gx/2) - 1) * Ax; end
if gy == 1,    starty = 2; else  starty = 2 + (floor(gy/2) - 1) * Ay; end

grapKernel(startx:startx+Ax-2, starty:starty+Ay-2) = 0.5;
grapKernel(startx:startx+Ax-1, starty:starty+Ay-2) = 0.5;
grapKernel(startx:startx+Ax-2, starty:starty+Ay-1) = 0.5;
%%%


% also declare variables for the coil-combined data
% (complex based on SVD weights, and magnitude from SSOS)
% This isn't needed for the final reconstruction, but can be really
% useful for debugging the motion-correction
if isempty(combinePars)   
    if isfield(twix_obj.hdr.MeasYaps,'asCoilSelectMeas') && strcmp(twix_obj.hdr.MeasYaps.asCoilSelectMeas{1}.asList{1}.sCoilElementID.tCoilID,'"32Ch_Head_7T"')
        % VE doesn't have 'asCoilSelectMeas', it's become sCoilSelectMeas and needs handling differently...
        disp('...............')
        disp('Detected use of 32Ch_Head_7T RF coil: therefore using predefined weights to also output a coil-combined image')
        disp('...............')
        combinePars = [-0.1049809 + -0.000000i;-0.01325038 - 0.1056198i;0.008624907 + 0.1050349i;-0.01387145 + 0.1499731i;0.009075958 - 0.2163172i;0.1473467 + 0.2160684i;0.2057536 - 0.003884733i;0.07486724 - 0.003189813i;-0.1640732 - 0.1206727i;-0.1645964 + 0.1266792i;-0.1423699 + 0.1314524i;-0.09954340 + 0.008023512i;0.2731563 + 0.1721494i;-0.2633240 - 0.2804754i;-0.2702133 + 0.05874889i;-0.07764274 + 0.05792535i;0.06793594 + 0.05160413i;-0.08903237 - 0.02729516i;0.0002752086 - 0.08664161i;-0.009266425 - 0.1078144i;-0.1005209 - 0.008859677i;0.1311229 - 0.05615813i;0.07896266 - 0.1288482i;0.02203858 - 0.06309662i;0.05616169 + 0.02123157i;-0.04540515 + 0.1785710i;-0.05066793 + 0.1244546i;0.03634416 + 0.003367372i;-0.1074157 + 0.1084992i;-0.2847729 + 0.04550047i;-0.1664367 + 0.06486133i;0.1344117 - 0.06306267i];
    elseif isfield(twix_obj.hdr.MeasYaps,'sCoilSelectMeas') && strcmp(twix_obj.hdr.MeasYaps.sCoilSelectMeas.aRxCoilSelectData{1}.asList{1}.sCoilElementID.tCoilID,'"Head_32"')
        disp('...............')
        disp('Detected use of Head_32 RF coil: therefore using predefined weights to also output a coil-combined image')
        disp('...............')
        combinePars = [0.1594946 + -0.000000i;0.1380453 - 0.1252680i;0.07394803 + 0.02236076i;-0.003868112 - 0.04217451i;0.1586883 + 0.07089107i;-0.01689743 - 0.1095780i;0.1480643 - 0.1408224i;0.1361057 - 0.2401873i;-0.1805168 - 0.007387850i;-0.03127036 - 0.06406011i;-0.01802163 - 0.06191686i;-0.1405896 + 0.07391697i;0.1038372 + 0.1059907i;-0.1225750 - 0.1220108i;0.1563336 + 0.02195730i;0.04546504 + 0.06252396i;-0.1738011 + 0.1473046i;-0.04153660 - 0.07626408i;-0.002118694 - 0.04707073i;-0.1151142 + 0.001793748i;0.09497871 + 0.1227228i;0.1143874 - 0.1130082i;0.1781285 + 0.1936812i;0.01422700 + 0.06057038i;-0.1706611 + 0.1959620i;-0.1682587 - 0.08187808i;-0.08768330 + 0.3116531i;-0.1621448 + 0.05920917i;-0.01707672 - 0.01533752i;-0.1027938 + 0.01049052i;0.3328240 - 0.1747643i;-0.1577251 - 0.09552535i];
    else
        combinePars = [];
    end
end

disp('...............')
disp('Declaring temporary variables to store GRAPPA recon')
tic

tempNameRoots.grappaRecon_1DFFT = [tempDir '/tempGRAPPAReconData_' MIDstr '_grappaRecon_1DFFT'];
tempNameRoots.reconSoS =          [tempDir '/tempGRAPPAReconData_' MIDstr '_reconSoS'];
tempFiles = dir([tempNameRoots.grappaRecon_1DFFT '*.mat']);
if ~isempty(tempFiles)
    delete([tempNameRoots.grappaRecon_1DFFT '*.mat']);
end
tempFiles = dir([tempNameRoots.reconSoS '*.mat']);
if ~isempty(tempFiles)
    delete([tempNameRoots.reconSoS '*.mat']);
end

if ~isempty(combinePars)
    tempNameRoots.dataCombined = [tempDir '/tempGRAPPAReconData_' MIDstr '_dataCombined'];
    tempFiles = dir([tempNameRoots.dataCombined '*.mat']);
    if ~isempty(tempFiles)
        delete([tempNameRoots.dataCombined '*.mat']);
    end
else
    tempNameRoots.dataCombined = [];
end



timingReport.declareVariables = toc; % this shouldn't take any time at all in the 'toDisk' version...!
disp('Done')
disp('...............')

    
%%
tic
disp('..............')
disp('...Now doing GRAPPA on each ''slice'' in the readout direction...')


parfor iS = 1:nread
    
    for iSet = 1:nSet
        
        for iEco = 1:nEco
            
            iPE2 = [1:nACS2]-floor(nACS2/2)+centerPar1;
            if any(iPE2<1 | iPE2>nACSmeas2) % seems this can happen for certain orientations / choices of parameters...
                iPE2 = 1:nACS2; % then hope that this works in this case...
            end
            
            this_fileid = fopen([tempDir '/temp_slicedata_' num2str(iS,'%.4d') '.bin'],'r');
            thisdata = fread(this_fileid,'single');
            fclose(this_fileid);
            if nSet > 1
                thisdata = reshape(thisdata,[nc npe1 nSet 2 npe2]);
                thisdata = squeeze(thisdata(:,:,iSet,:,:));
            elseif nEco > 1
                thisdata = reshape(thisdata,[nc npe1 nEco 2 npe2]);
                thisdata = squeeze(thisdata(:,:,iEco,:,:));
            else
                thisdata = reshape(thisdata,[nc npe1 2 npe2]);
            end
            thisdata = squeeze(thisdata(:,:,1,:) + 1i*thisdata(:,:,2,:));
            thisdata = permute(thisdata,[2 3 1]);
            
            startX = find(thisdata(:,1,1),1,'first');
            
            if nSliceNeighbours > 0
                
                src = []; targ = [];
                
                for this_iS = iS-nSliceNeighbours:iS+nSliceNeighbours
                    if this_iS >= 1 && this_iS <= nread
                        thisACSdataFull = squeeze(ACSdata(this_iS,:,:,:,1,1,1,iEco,1,iSet));
                        thisACSdataFull = permute(thisACSdataFull,[2 3 1]);
                        thisACSdata = thisACSdataFull([1:nACS1]-floor(nACS1/2)+floor(nACSmeas1/2),iPE2,:);
                        [~, this_src, this_targ] = GrappaCalib3D_arb(thisACSdata,grapKernel,0,0);
                        src = [src; this_src];
                        targ = [targ; this_targ];
                    end
                end
                
                % and store the current slice for reinsertion
                thisACSdataFull = squeeze(ACSdata(iS,:,:,:,1,1,1,iEco,1,iSet));
                thisACSdataFull = permute(thisACSdataFull,[2 3 1]);
                
                thisGrapW = pinv(src)*targ;
                
            else
                
                thisACSdataFull = squeeze(ACSdata(iS,:,:,:,1,1,1,iEco,1,iSet));
                thisACSdataFull = permute(thisACSdataFull,[2 3 1]); % this permute is necessary for my GRAPPA code, but may not be the fastest approach...
                thisACSdata = thisACSdataFull([1:nACS1]-floor(nACS1/2)+floor(nACSmeas1/2),iPE2,:);
                
                thisGrapW = GrappaCalib3D_arb(thisACSdata,grapKernel,0,1);
            end
            
            thisGrapRecon = GrappaReco3D_arb(thisdata,grapKernel,thisGrapW,[Ax 1],[startX 1]);
            
            
            %%
            % and put the ACS lines back in
            % WARNING - I don't know if the following line will be correct for all
            % possible image resolution/matrix sizes, so possibly safer to leave this to the
            % separate image recon part...
            thisGrapRecon([1:nACSmeas1]-floor(nACSmeas1/2)+centerLin1-1,:,:) = thisACSdataFull;
            
            if ~isempty(combinePars)
                thisCombination = complex(zeros(npe1,npe2,'single'));
                for iC = 1:nc
                    thisCombination = thisCombination + combinePars(iC)*thisGrapRecon(:,:,iC);
                end
%                 dataCombined(iS,:,:,iSet) = reshape(thisCombination,[1 npe1 npe2]);
                
                parsave([tempNameRoots.dataCombined '_' num2str(iS) '_' num2str(iEco) '_' num2str(iSet) '.mat'],thisCombination);
            end
            
            thisImage = ssos(ifft2s(thisGrapRecon));
            %                 mOutGRAPPA.reconSoS(iS,:,:,iSet) = reshape(thisImage,[1 npe1 npe2]);
            %                 mOutGRAPPA.grappaRecon_1DFFT(iS,:,:,:,iSet) = reshape(permute(thisGrapRecon,[3 1 2]),[1 nc npe1 npe2]);
            outData = reshape(thisImage,[1 npe1 npe2]);
            parsave([tempNameRoots.reconSoS '_' num2str(iS) '_' num2str(iEco) '_' num2str(iSet) '.mat'],outData);
            outData = reshape(permute(thisGrapRecon,[3 1 2]),[1 nc npe1 npe2]);
            parsave([tempNameRoots.grappaRecon_1DFFT '_' num2str(iS) '_' num2str(iEco) '_' num2str(iSet) '.mat'],outData);
            
        end % iEco loop
        
    end % iSet loop        
    
end % iS for 'slices' in readout direction



%% recompile the reconSoS and dataCombined into a single file as these should be small in total
reconSoS = zeros(nread,npe1,npe2,nSet,nEco,'single');
if ~isempty(combinePars)
    dataCombined = complex(zeros(nread,npe1,npe2,nSet,nEco,'single'));
end
for iS = 1:nread   
    for iSet = 1:nSet        
        for iEco = 1:nEco
            
            thisData = load([tempNameRoots.reconSoS '_' num2str(iS) '_' num2str(iEco) '_' num2str(iSet) '.mat']);
            reconSoS(iS,:,:,iSet,iEco) = thisData.outData;
            
           if ~isempty(combinePars)
                thisData = load([tempNameRoots.dataCombined '_' num2str(iS) '_' num2str(iEco) '_' num2str(iSet) '.mat']);
                dataCombined(iS,:,:,iSet,iEco) = thisData.thisCombination;
            end
            
        end
    end
end
delete([tempNameRoots.reconSoS '*.mat']);
save(tempNameRoots.reconSoS,'reconSoS','-v7.3');

if ~isempty(combinePars)
    delete([tempNameRoots.dataCombined '*.mat']);
    save(tempNameRoots.dataCombined,'dataCombined','-v7.3');
end


fprintf('\n');
disp('Done')
disp('..............')
timingReport.GRAPPArecon = toc;


%% delete temporary files which are currently not needed...

delete([tempDir '/temp_slicedata_*.bin'])
