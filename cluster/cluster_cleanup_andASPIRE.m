load(cleanupFile)

run([RETROMOCOBOX_PATH '/addRetroMoCoBoxToPath.m']);
addpath(SPM_PATH);

run([ASPIRE_HOME '/aspire_startup.m']);

%% Run ASPIRE on both noMoCo and MoCo




data.read_dir = tempDir;
data.filename_mag = fullfile(data.read_dir, 'GRE_5D_mag.nii'); % 5D NIfTI input (x,y,z,eco,cha)
data.filename_phase = fullfile(data.read_dir, 'GRE_5D_phs.nii');
data.write_dir = fullfile(tempDir,'ASPIRE_noMoCo');

% OPTIONS
data.poCalculator = AspirePoCalculator; % AspireBipolarPoCalculator('non-linear correction') for bipolar acquisitions (at least 3 echoes)
data.smoother = NanGaussianSmoother; % NanGaussianSmoother, GaussianBoxSmoother (=default)

% data.processing_option = 'all_at_once'; % all_at_once, slice_by_slice (slice_by_slice requires fslmerge)
data.processing_option = 'slice_by_slice';

% run ASPIRE
tic;
run(Aspire(data));
aspireTime = secs2hms(toc); 

disp(['Whole calculation took: ' aspireTime]);
disp(['Files written to: ' data.write_dir]);

data2 = data;
data2.filename_mag = fullfile(data2.read_dir,'GRE_5D_MoCo_mag.nii');
data2.filename_phase = fullfile(data2.read_dir,'GRE_5D_MoCo_phs.nii');
data2.write_dir = fullfile(tempDir,'ASPIRE_MoCo');

% run ASPIRE
tic;
run(Aspire(data2));
aspireTime = secs2hms(toc); 

disp(['Whole calculation took: ' aspireTime]);
disp(['Files written to: ' data2.write_dir]);

%% Copy ASPIRE output to be main GRE recon

copyfile(fullfile(data.write_dir,'results','combined_mag.nii'),fullfile(outDir,'GRE_mag.nii'));
copyfile(fullfile(data.write_dir,'results','combined_phase.nii'),fullfile(outDir,'GRE_phs.nii'));
copyfile(fullfile(data2.write_dir,'results','combined_mag.nii'),fullfile(outDir,'GRE_MoCo_mag.nii'));
copyfile(fullfile(data2.write_dir,'results','combined_phase.nii'),fullfile(outDir,'GRE_MoCo_phs.nii'));


%%

fidHTML = fopen([htmlDir '/index.html'],'a');



%% And put the reconstructed images into the html

all_ims = rn([outDir '/GRE_mag.nii']);
all_ims = all_ims(:,:,:,1); % keep just first echo
all_ims_corrected = rn([outDir '/GRE_MoCo_mag.nii']);
all_ims_corrected = all_ims_corrected(:,:,:,1);

ov1 = orthoview(all_ims,'drawIms',0);
imab_overwrite([htmlDir '/GRE_mag.png'],ov1.oneIm);
ov1 = orthoview(all_ims_corrected,'drawIms',0);
imab_overwrite([htmlDir '/GRE_MoCo_mag.png'],ov1.oneIm);

fprintf(fidHTML,['GRE image before correction:<br>\n']);
fprintf(fidHTML,['<img src="GRE_mag.png"><br><br>\n']);
fprintf(fidHTML,['GRE image after correction:<br>\n']);
fprintf(fidHTML,['<img src="GRE_MoCo_mag.png"><br><br>\n']);


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
% fprintf(fidHTML,['<strong>Application of retrospective motion-correction: </strong>' num2str(round(timingReport_totalTimeApplyMoco/60)) ' minutes\n']); % <--- not passed through currently in cluster mode

% include version number
fprintf(fidHTML,['<br><br><br><em>' char(datetime) '- created with reconstructSiemensGREwithFatNavs.m, version: ' reconPars.retroMocoBoxVersion '</em>\n']);


fprintf(fidHTML,'</body></html>\n');
fclose(fidHTML);


%% Delete the temporary files (which could be rather large...!)
 
% clear temporary
if ~reconPars.bKeepGRAPPArecon
    fclose('all'); % ASPIRE currently doesn't close files properly
    delete([tempDir '/GRE_*.nii'])
    rmdir([tempDir '/matlab_tmp*'],'s')
    rmdir([tempDir '/ASPIRE_*'],'s');
    
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
