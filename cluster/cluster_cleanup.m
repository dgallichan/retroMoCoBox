load(cleanupFile)

run([RETROMOCOBOX_PATH '/addRetroMoCoBoxToPath.m']);
run([MIRT_PATH '/setup.m']);
addpath(SPM_PATH);



%%

fidHTML = fopen([htmlDir '/index.html'],'a');


tFinish_applyMoco = clock;
timingReport_totalTimeApplyMoco = etime(tFinish_applyMoco,tStart_applyMoco);

%% And put the reconstructed images into the html

all_ims = rn([outDir '/GRE.nii']);
all_ims_corrected = rn([outDir '/GRE_corrected.nii']);

ov1 = orthoview(all_ims,'drawIms',0);
imab_overwrite([htmlDir '/GRE.png'],ov1.oneIm);
ov1 = orthoview(all_ims_corrected,'drawIms',0);
imab_overwrite([htmlDir '/GRE_corrected.png'],ov1.oneIm);

fprintf(fidHTML,['GRE image before correction:<br>\n']);
fprintf(fidHTML,['<img src="GRE.png"><br><br>\n']);
fprintf(fidHTML,['GRE image after correction:<br>\n']);
fprintf(fidHTML,['<img src="GRE_corrected.png"><br><br>\n']);


testMagick = system('convert -version');

if testMagick==0 % can use ImageMagick to make animated GIFs...
    processString = ['convert -dispose 2 -delay 50 -loop 0 ' htmlDir '/GRE*.png ' htmlDir '/mov_GRE.gif'];
    system(processString);
    fprintf(fidHTML,['GRE image movie before/after correction:<br>\n']);
    fprintf(fidHTML,['<img src="mov_GRE.gif"><br><br>\n']);
end



%% And put a timing report into the html

stopTime = clock;
totalTime = etime(stopTime,startTime)/60/60;
totalTime_hrs = floor(totalTime);
if totalTime_hrs > 0
    totalTime_mins = round(rem(totalTime,totalTime_hrs)*60);
else
    totalTime_mins = round(totalTime*60);
end

fprintf(fidHTML,['<h4>Total reconstruction time: ' num2str(totalTime_hrs) ' hours, ' num2str(totalTime_mins) ' mins</h4>\n']);
fprintf(fidHTML,['<strong>Calculate GRAPPA weights for FatNavs: </strong>' num2str(round(timingReport_FatNavs.calculateGRAPPAweights)) ' seconds.<br>\n']);
fprintf(fidHTML,['<strong>Reconstruct FatNavs: </strong>' num2str(nFatNavs) 'x ' num2str(round(mean(timingReport_FatNavs.eachFatNav))) ' seconds. Total time (possibly parallelized!): = ' num2str(round(timingReport_FatNavs.allFatNavs)) ' seconds. <br>\n']);
fprintf(fidHTML,['<strong>Align FatNavs using SPM: </strong>' num2str(round(timingReport_FatNavs.SPMalignment)) ' seconds.<br>\n']);
fprintf(fidHTML,['<strong>Application of retrospective motion-correction: </strong>' num2str(round(timingReport_totalTimeApplyMoco/60)) ' minutes\n']);

% include version number
fprintf(fidHTML,['<br><br><br><em>' char(datetime) '- created with reconstructSiemensGREwithFatNavs.m, version: ' retroMocoBoxVersion '</em>\n']);


fprintf(fidHTML,'</body></html>\n');
fclose(fidHTML);


%% Delete the temporary files (which could be rather large...!)
 
% clear temporary
if ~reconPars.bKeepGRAPPArecon
    if ~isempty(tempNameRoots.grappaRecon_1DFFT)
        delete([tempNameRoots.grappaRecon_1DFFT '*.mat']);
    end
    if ~isempty(tempNameRoots.reconSoS)
        delete([tempNameRoots.reconSoS '*.mat']);
    end
    if ~isempty(tempNameRoots.dataCombined)
        delete([tempNameRoots.dataCombined '*.mat']);
    end
    if ~isempty(tempNameRoots.allReconPars)
        delete([tempNameRoots.allReconPars]);
    end
    if ~isempty(tempNameRoots.finalImage)
        delete([tempNameRoots.finalImage '*.mat']);                
    end
    if ~isempty(tempNameRoots.cleanupFiles)
        delete([tempNameRoots.cleanupFiles]);
    end
    if ~isempty(tempNameRoots.clusterScript)
        delete([tempNameRoots.clusterScript]);
    end    
    if ~isempty(tempNameRoots.clusterScriptRecombine)
        delete([tempNameRoots.clusterScriptRecombine]);
    end
    if ~isempty(tempNameRoots.clusterScriptCleanup)
        delete([tempNameRoots.clusterScriptCleanup]);
    end
    try
        rmdir(tempDir)
    catch
        disp(['Warning: Couldn''t delete temporary folder "' tempDir '" - folder not empty?'])
    end
end

if reconPars.bZipNIFTIs
    gzip([outDir '/*.nii']);
    delete([outDir '/*.nii']);
end

if ~reconPars.bKeepFatNavs
    rmdir(fatnavdir,'s')
end
    


fprintf('*************************************************************\n')
fprintf('***** reconstructSiemensGREwithFatNavs.m completed! *****\n')
fprintf('*************************************************************\n')
fprintf(['Total reconstruction time: ' num2str(totalTime_hrs) ' hours, ' num2str(totalTime_mins) ' mins\n']);
