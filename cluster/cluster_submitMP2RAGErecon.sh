#!/bin/bash

if [ $# -lt 1 ];then
	echo ""
	echo "Runs MP2RAGE with FatNavs recon on the CUBRIC cluster."
	echo ""
	echo "Usage: cluster_submitMP2RAGErecon.sh <input_raw_twix_file>"
	echo ""
	echo "IMPORTANT:"
	echo "Before running you will need to make sure that you set"
	echo "the paths to the retroMoCoBox and SPM"
	echo "correctly inside this script."
	echo ""
	exit 1;
fi

inputfile=$1

export RETROMOCOBOX_HOME=/home/scedg10/retroMoCoBox
export SPM_HOME=/cubric/software/spm.versions/spm12
export CLUSTER_LOG_PATH=${PWD}/clusterlogs

# if 'clusterlogs' folder doesn't exist already, create it:
[ -d clusterlogs ] || mkdir clusterlogs

tempfile="$(mktemp)"
echo "#!/bin/bash" > ${tempfile}
echo "#SBATCH --job-name=MP2RAGEFatNavs" >> ${tempfile}
echo "#SBATCH -p cubric-default" >> ${tempfile}
echo "#SBATCH -o ${CLUSTER_LOG_PATH}/MP2RAGE_FatNavs_%j.out" >> ${tempfile}
echo "#SBATCH -e ${CLUSTER_LOG_PATH}/MP2RAGE_FatNavs_%j.err" >> ${tempfile}
echo "#SBATCH --mem=160G" >> ${tempfile}
echo "#SBATCH --cpus-per-task 10" >> ${tempfile}
echo "#SBATCH --ntasks 1" >> ${tempfile}
echo "cd ${RETROMOCOBOX_HOME}" >> ${tempfile}
echo "matlab -nodisplay -nodesktop -nosplash -r \"rawDataFile='${PWD}/${inputfile}';swapDims_xyz = [1 0 1]; FatNavRes_mm = 2;run_SiemensMP2RAGErecon_cluster;exit;\"" >> ${tempfile}

echo " "
echo " "
echo "Creating this temporary file to process your data:"
echo "=================================================="
cat ${tempfile}
echo " "
echo " "

sbatch ${tempfile}
rm ${tempfile}

