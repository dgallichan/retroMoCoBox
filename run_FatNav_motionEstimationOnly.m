% Process the FatNavs and export motion parameters in the coordinate frame
% of the host sequence acquisition - useful for using alternative to
% retroMoCoBox for the MoCo processing

run('~/retroMoCoBox/addRetroMoCoBoxToPath.m')

rawDataFile = '/cubric/data/scedg10/exampleFatNavsRawData/meas_MID101_mp2rage_FatNav_1mm_smallMotion.dat';
outRoot = '/home/scedg10/myscratch/retroMoCoBox_unitTests/tests/';


%%

twix_obj = mapVBVD_fatnavs(rawDataFile,'removeOS',1);

processFatNavs_GRAPPA4x4(twix_obj, outRoot, 'bSwapHandedness',1);


%%

rotAndShift = getSiemensRotMatAndShift(twix_obj.hdr);

MIDstr = getMIDstr(rawDataFile);
fitResult = load([outRoot '/motion_parameters_spm_' MIDstr '.mat']);
    
% conver the mpars into coordinate frame of host sequence
% (corresponds to default code from retroMoCoBox v1.0.2)
this_fitMat = fitResult.MPos_cent.mats;

extraFlipMat1 = diag([-1 1 -1]);
extraFlipMat2 = eye(4);
extraPositionOffsetSignFlips = [1 -1 1];

A_fatnav2host = eye(4);
% add rotations:
A_fatnav2host(1:3,1:3) = rotAndShift.RotMat*extraFlipMat1;
% and translations:
A_fatnav2host(1:3,4) = -rotAndShift.RotMat*extraFlipMat1*...
    (extraPositionOffsetSignFlips(:).*rotAndShift.Shifts_SagCorTra(:));

A_fatnav2host_forMats = extraFlipMat2*A_fatnav2host;

A_mpars_mm = moveFrame(this_fitMat,A_fatnav2host_forMats);

%%

save([outRoot '/motion_parameters_spm_' MIDstr '_inHostFOV.mat'],'A_mpars_mm')

