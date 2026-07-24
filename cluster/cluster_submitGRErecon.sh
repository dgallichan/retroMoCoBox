#!/bin/bash

if [ $# -lt 1 ];then
	echo ""
	echo "Runs GRE with FatNavs recon on the CUBRIC cluster."
	echo ""
	echo "Usage: cluster_submitGRErecon.sh <full_path_to_input_raw_twix_file> <full_path_to_output_root_folder>"
	echo ""
	echo "If you pass just one input file you will get the outputs in the current working folder."
	echo ""
	echo "A folder called 'GRErecon_MIDXXX' will be created in the output root folder,"
	echo "where the XXX correponds to the MID number of the raw data file. A folder called"
	echo "'clusterlogs' will also be created (if it doesn't already exist) which stores"
	echo "the log files from the jobs running on the cluster (useful for debugging if"
	echo "something should go wrong!)"
	echo "While the job is running, there will also be a couple of other temporary"
	echo "folders created - which can run into >200 GB of data usage for high-res"
	echo "multi-echo data. Please make sure you have enough quota to handle this,"
    echo "especially if attempting to run multiple reconstructions simultaneously!"
	echo ""
	echo "IMPORTANT:"
	echo "Before running you will need to make sure that you set"
	echo "the RETROMOCOBOX_HOME environment variable and also ASPIRE_HOME"
	echo "by manually editing this script and"
	echo "specifying them below this text in the .sh file itself."
	exit 1;
fi

export RETROMOCOBOX_HOME=/cubric/software/matlab.toolboxes/retroMoCoBox_v1.0.1
# export RETROMOCOBOX_HOME=/home/scedg10/retroMoCoBox
export ASPIRE_HOME=/cubric/software/matlab.toolboxes/ASPIRE
export SPM_HOME=/cubric/software/spm.versions/spm12

# The GRE recon is more complicated that the MP2RAGE because at high resolutions with several echoes
# the RAM requirements for the 'all in RAM' option are quickly overwhelming.
#
# This script therefore initiates a series of scripts to chunk the problem to reduce RAM requirements
# but this comes with the drawback of needing lots of disk space (often >200 GB) for temporary files
# that get shared across separate cluster jobs.
#
# The current pathway for recon that this takes:
#
# Initial cluster job (job name: "GRE_FatNavs")
#     1. run_SiemensGRErecon_cluster.m 
#          - this is just a helper script that checks paramters, then launches reconstructSiemensGREwithFatNavs_cluster.m
#     2. reconstructSiemensGREwithFatNavs_cluster.m
#          - this processes the FatNavs and does GRAPPA on the host GRE sequence, then spawns 3 further cluster jobs for the rest
#
# Array cluster job (job name: "GREreconHelperArray" - one job per echo):
#     cluster_runMultiEchoGRE_eachEcho.m 
#         - this is run on each echo of the GRE separately to apply the MoCo
#
# Coil combination job (job name: "GREreconRecombine"):
#     cluster_combineCoils_forASPIRE.m  
#         - this makes (large!) 5D NIFTI files with separate magnitude and phase in the format that ASPIRE code expects
#
# Cleanup and post-processing (ASPIRE) job (job name: "GREreconCleanup"):
#     cluster_cleanup_andASPIRE.m 
#         - this runs ASPIRE (if more than one TE, otherwise sum-of-squares)
#         - then finishes off HTML report and deletes temporary files and folder
#
# For 0.67mm data with 7 echoes (20GB raw data file), this takes ~4hrs to run on the cluster and requires 209 GB of additional disk 
# space while running.
#
# Note that if you still have a sub-folder starting 'temp_' in the output folder then either the reconstruction has not yet finished
# or there has been an error with the processing. Check if you still have jobs queued on the cluster, and check the 'clusterlogs' 
# folder for error messages.


inputfile=$1

if [ $# -lt 2 ];then
	outRoot=${PWD}
    export CLUSTER_LOG_PATH=${PWD}/clusterlogs
else
	outRoot=$2
    export CLUSTER_LOG_PATH=${outRoot}/clusterlogs
fi



# if 'cluster log path' folder doesn't exist already, create it:
[ -d ${CLUSTER_LOG_PATH} ] || mkdir ${CLUSTER_LOG_PATH}

tempfile="$(mktemp)"

echo "#!/bin/bash" > ${tempfile}
echo "#SBATCH --job-name=GREFatNavs" >> ${tempfile}
# echo "#SBATCH -p cubric-centos7" >> ${tempfile}
echo "#SBATCH -p cubric-rocky8" >> ${tempfile}
echo "#SBATCH -o ${CLUSTER_LOG_PATH}/GRE_FatNavs_%j.out" >> ${tempfile}
echo "#SBATCH -e ${CLUSTER_LOG_PATH}/GRE_FatNavs_%j.err" >> ${tempfile}
echo "#SBATCH --mem=80G" >> ${tempfile}
echo "#SBATCH --cpus-per-task 12" >> ${tempfile}
echo "#SBATCH --ntasks 1" >> ${tempfile}
echo "cd ${RETROMOCOBOX_HOME}" >> ${tempfile}
echo "matlab -nodisplay -nodesktop -nosplash -r \"rawDataFile='${inputfile}';outRoot='${outRoot}';bKeepGRAPPArecon=0;parpoolSize = 6;script_pwd = '${PWD}';ASPIRE_HOME='${ASPIRE_HOME}';run_SiemensGRErecon_cluster;exit;\"" >> ${tempfile}
echo " "
echo " "
echo "Creating this temporary file to process your data:"
echo "=================================================="
echo " "
cat ${tempfile}
echo " "
echo " "

sbatch ${tempfile}

rm ${tempfile}
