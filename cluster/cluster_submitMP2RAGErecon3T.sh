#!/bin/bash

if [ $# -lt 1 ];then
	echo ""
	echo "Runs MP2RAGE with FatNavs recon on the CUBRIC cluster with parameters for 3T acquisitions."
	echo ""
	echo "Usage: cluster_submitMP2RAGErecon3T.sh <full_path_to_input_raw_twix_file> <full_path_to_output_root_folder>"
	echo ""
	echo "If you pass just one input file you will get the outputs in the current working folder."
	echo ""
	echo ""
	echo "IMPORTANT:"
	echo "Before running you will need to make sure that you set"
	echo "the RETROMOCOBOX_HOME environment variable manually within this script itself!"
	echo ""
	exit 1;
fi

# Note that at the moment the only thing that really changes between 7T and 3T is that the FatNav resolution will be 4mm instead of 2mm
# This is set manually in the line that calls MATLAB below.

inputfile=$1

export RETROMOCOBOX_HOME=/home/scedg10/retroMoCoBox

export SPM_HOME=/cubric/software/spm.versions/spm12

if [ $# -lt 2 ];then
	outRoot=${PWD}
    export CLUSTER_LOG_PATH=${PWD}/clusterlogs
else
	outRoot=$2
    export CLUSTER_LOG_PATH=${outRoot}/clusterlogs
fi

# if 'clusterlogs' folder doesn't exist already, create it:
[ -d clusterlogs ] || mkdir clusterlogs

tempfile="$(mktemp)"
echo "#!/bin/bash" > ${tempfile}
echo "#SBATCH --job-name=MP2RAGEFatNavs" >> ${tempfile}
echo "#SBATCH -p cubric-rocky8" >> ${tempfile}
echo "#SBATCH -o ${CLUSTER_LOG_PATH}/MP2RAGE_FatNavs_%j.out" >> ${tempfile}
echo "#SBATCH -e ${CLUSTER_LOG_PATH}/MP2RAGE_FatNavs_%j.err" >> ${tempfile}
echo "#SBATCH --mem=160G" >> ${tempfile}
echo "#SBATCH --cpus-per-task 10" >> ${tempfile}
echo "#SBATCH --ntasks 1" >> ${tempfile}
echo "cd ${RETROMOCOBOX_HOME}" >> ${tempfile}
echo "matlab -nodisplay -nodesktop -nosplash -r \"rawDataFile='${PWD}/${inputfile}';swapDims_xyz = [0 0 1]; FatNavRes_mm = 4;run_SiemensMP2RAGErecon_cluster;exit;\"" >> ${tempfile}

echo " "
echo " "
echo "Creating this temporary file to process your data:"
echo "=================================================="
cat ${tempfile}
echo " "
echo " "


sbatch ${tempfile}
rm ${tempfile}

