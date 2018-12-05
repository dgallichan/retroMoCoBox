function [grappaRecon_1DFFT, mOutGRAPPA] = performHostGRAPPArecon_RAMonly(twix_obj, counterStruct, nSliceNeighbours, combinePars, lambda)
%
% function [grappaRecon_1DFFT, mOutGRAPPA] = performHostGRAPPArecon_RAMonly(twix_obj, counterStruct, nSliceNeighbours, combinePars, lambda)
%
% Applies GRAPPA to the host data, which here is assumed to have acceleration in
% the first phase-encoding direction only (as is the current limitation of
% the MP2RAGE Siemens code).
%
% This version attempts to loop over 'Echoes' and 'Sets'. 
%
% If you want to use other counters (currently only 'Repetitions' or 'Averages') then 
% specify which to reconstruct with counterStruct.iRep or
% counterStruct.iAve. It is assumed that there is separate ACS data available for each repetition.
%
%
% -------
% daniel.gallichan@epfl.ch, June 2015
%
% 11/3/16 - danielg - created '_RAMonly' version which should handle more
%                     arbitrary datasets, as when using matfiles it is
%                     difficult to handle arrays which might have a different
%                     number of dimensions. Also removed the coil
%                     compression as it doesn't seem to be compatible with
%                     motion-correction. As it's RAM only the mOutGRAPPA
%                     structure is also only constructed at the end, as then 
%                     it becomes compatible with parfor...
%
% 5/2/18 - gallichand@cardiff.ac.uk
%        - Introduce the 'counterStruct' to attempt to handle repetitions
%          and averages as simply as possible. NB: Should now handle
%          'averages' properly - but not 'repetitions' as the counters in
%          the FatNav recon will then clash. This is fixable, but not yet
%          implemented...
   
if nargin < 5
    lambda = 0; % regularisation term for determination of GRAPPA weights
end

if nargin < 4
    combinePars = [];
end

if nargin < 3
    nSliceNeighbours = 2; % no. of neighbouring 'virtual' slices to include in ACS lines to try to improve conditioning and regularize in the readout direction
    % quick tests suggest that 2 is a good number to use... :)
end

if nargin < 2
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


%%%%%%%%%%%%%%%%%%%%
%%% Will only currently work with 1D GRAPPA...!
%%%%%%%%%%%%%%%%%%%%
%%


figIndex = 999; % use this figure number for live 'show' of how things are going...

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

tic
disp('..............')
disp('... Loading the raw data file into RAM and doing a 1D-FFT')
allData = twix_obj.image(:,:,:,:,:,counterStruct.iAve,:,:,counterStruct.iRep,:);
timingReport.FFTperSlice = toc;
allData = ifft1s(allData,1);

disp('Done')
disp('..............')


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
grappaRecon_1DFFT = complex(zeros(nread,nc,npe1,npe2,nEco,nSet,'single'));
reconSoS = zeros(nread,npe1,npe2,nEco,nSet,'single');
if ~isempty(combinePars)
    dataCombined = complex(zeros(nread,npe1,npe2,nEco,nSet,'single'));
end
timingReport.declareVariables = toc;
disp('Done')
disp('...............')

    
%%
tic
disp('..............')
disp('...Now doing GRAPPA on each ''slice'' in the readout direction...')

parfor iS = 1:nread 
    
    for iEco = 1:nEco
        
        for iSet = 1:nSet
            
            
            
            iPE2 = [1:nACS2]-floor(nACS2/2)+centerPar1;
            if any(iPE2<1 | iPE2>nACSmeas2) % seems this can happen for certain orientations / choices of parameters...
                iPE2 = 1:nACS2; % then hope that this works in this case...
            end
                        
            thisdata = squeeze(allData(iS,:,:,:,1,1,1,iEco,1,iSet));           
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
                
                if (lambda)
                    lambdaApply = lambda*norm(src(:));
                    thisGrapW  = pinv(src'*src + lambdaApply^2*eye(size(src,2))) * src' * targ;
                else
                    thisGrapW = pinv(src)*targ;
                end
                          
            else
                thisACSdataFull = squeeze(ACSdata(iS,:,:,:,1,1,1,iEco,1,iSet));
                thisACSdataFull = permute(thisACSdataFull,[2 3 1]); % this permute is necessary for my GRAPPA code, but may not be the fastest approach...
                thisACSdata = thisACSdataFull([1:nACS1]-floor(nACS1/2)+floor(nACSmeas1/2),iPE2,:);
                
                thisGrapW = GrappaCalib3D_arb(thisACSdata,grapKernel,lambda,1);
            end           
            
            
            thisGrapRecon = GrappaReco3D_arb(thisdata,grapKernel,thisGrapW,[Ax 1],[startX 1]);
            
            
            %
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
                dataCombined(iS,:,:,iEco,iSet) = reshape(thisCombination,[1 npe1 npe2]);
             end
            
            thisImage = ssos(ifft2s(thisGrapRecon));
            
            
   %%
            reconSoS(iS,:,:,iEco,iSet) = reshape(thisImage,[1 npe1 npe2]);
            grappaRecon_1DFFT(iS,:,:,:,iEco,iSet) = reshape(permute(thisGrapRecon,[3 1 2]),[1 nc npe1 npe2]);

            
%             if ~isempty(figIndex)
%                 fig(figIndex)
%                 clf
%                 imab(thisImage)
%                 title(['GRAPPA recon, slice ' num2str(iS) ' out of ' num2str(nread) ', echo index: ' num2str(iEco) ', set index: ' num2str(iSet)])
%                 drawnow
%             end
            
            
        end % iSet loop
        
    end
    
    if mod(iS,10)==0
        fprintf('.');
    end
    
end % iS for 'slices' in readout direction
fprintf('\n');
disp('Done')
disp('..............')
timingReport.GRAPPArecon = toc;

mOutGRAPPA.reconSoS = reconSoS;
mOutGRAPPA.timingReport = timingReport;
if ~isempty(combinePars)
    mOutGRAPPA.dataCombined = dataCombined;
end


