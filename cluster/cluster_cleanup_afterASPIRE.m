load(cleanupFile)

run([RETROMOCOBOX_PATH '/addRetroMoCoBoxToPath.m']);
% run([MIRT_PATH '/setup.m']);
addpath(SPM_PATH);



%%

fidHTML = fopen([htmlDir '/index.html'],'a');



%% And put the reconstructed images into the html

all_ims = rn([outDir '/GRE_mag']);
all_ims = all_ims(:,:,:,1); % keep just first echo
all_ims_corrected = rn([outDir '/GRE_MoCo_mag']);
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

Hxyz = size(all_ims);

%% Add a zoom of the GRE images before and after correction
xi1 = round(Hxyz(1)*.5:Hxyz(1)*.75);
yi1 = round(Hxyz(2)*.7: .9*Hxyz(2)); % arbitrary cut off points for zoom in x,y and z
zi1 = round(Hxyz(3)/2);

xi2 = round(Hxyz(1)*.55:.8*Hxyz(1));
yi2 = round(Hxyz(2)/2); % arbitrary cut off points for zoom in x,y and z
zi2 = round(Hxyz(3)*.6: .95*Hxyz(3));

xi3 = round(Hxyz(1)/2+Hxyz(1)/8);
yi3 = round(Hxyz(2)/2+Hxyz(2)/8: .9*Hxyz(2)); % arbitrary cut off points for zoom in x,y and z
zi3 = round(Hxyz(3)/2: .9*Hxyz(3));



clim1 = percentile(all_ims,99);
clim2 = percentile(all_ims_corrected,99);

figure
set(gcf,'Position',[           122         592        1665         517])
subplot1(1,3)
subplot1(1)
imab(squeeze(all_ims(xi1,yi1,zi1)),[0 clim1])
subplot1(2)
imab(squeeze(all_ims(xi2,yi2,zi2)),[0 clim1])
subplot1(3)
imab(squeeze(all_ims(xi3,yi3,zi3)),[0 clim1])
colormap(gray)

export_fig([htmlDir '/zoom.png'])
clf
subplot1(1,3)
subplot1(1)
imab(squeeze(all_ims_corrected(xi1,yi1,zi1)),[0 clim2])
subplot1(2)
imab(squeeze(all_ims_corrected(xi2,yi2,zi2)),[0 clim2])
subplot1(3)
imab(squeeze(all_ims_corrected(xi3,yi3,zi3)),[0 clim2])
colormap(gray)
export_fig([htmlDir '/zoom_corrected.png'])

if testMagick==0
    processString = ['convert -dispose 2 -delay 50 -loop 0 ' htmlDir '/zoom.png ' htmlDir '/zoom_corrected.png ' htmlDir '/mov_zoom.gif'];
    system(processString);
    
    fprintf(fidHTML,['Zooms before/after correction:<br>\n']);
    fprintf(fidHTML,['<img src="mov_zoom.gif"><br><br>\n']);
else
    fprintf(fidHTML,['Zooms before correction:<br>\n']);
    fprintf(fidHTML,['<img src="zoom.png"><br><br>\n']);
    fprintf(fidHTML,['Zooms after correction:<br>\n']);
    fprintf(fidHTML,['<img src="zoom_corrected.png"><br><br>\n']);
end

%% Add zooms of minIPs of centre

xi = round(.3*Hxyz(1):.7*Hxyz(1));
yi = round(.3*Hxyz(2):.7*Hxyz(2));
zi = round(.4*Hxyz(3):.7*Hxyz(3));

% figure
% imab(min(all_ims(xi,yi,zi),[],3))
% colormap(gray)

imab_overwrite([htmlDir '/GRE_zoom_minIP.png'],imresize(min(all_ims(xi,yi,zi),[],3),3));
imab_overwrite([htmlDir '/GRE_zoom_corrected_minIP.png'],imresize(min(all_ims_corrected(xi,yi,zi),[],3),3));

if testMagick==0
    processString = ['convert -dispose 2 -delay 50 -loop 0 ' htmlDir '/GRE_zoom_minIP.png ' htmlDir '/GRE_zoom_corrected_minIP.png ' htmlDir '/mov_zoom_minIP.gif'];
    system(processString);
    
    fprintf(fidHTML,['Zooms of minIP before/after correction:<br>\n']);
    fprintf(fidHTML,['<img src="mov_zoom_minIP.gif"><br><br>\n']);
else
    fprintf(fidHTML,['Zooms of minIP before correction:<br>\n']);
    fprintf(fidHTML,['<img src="GRE_zoom_minIP.png"><br><br>\n']);
    fprintf(fidHTML,['Zooms of minIP after correction:<br>\n']);
    fprintf(fidHTML,['<img src="GRE_zoom_corrected_minIP.png"><br><br>\n']);
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
    
    delete([tempDir '/GRE_*.nii']);
    rmdir([tempDir '/combined_MoCo'],'s')
    rmdir([tempDir '/combined_noMoCo'],'s')

    
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

if reconPars.bZipNIFTIs && ~isempty(dir([outDir '/*.nii']))
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
