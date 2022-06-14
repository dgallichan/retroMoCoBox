#!/bin/bash

if [ $# -lt 1 ];then
	echo ""
	echo "Runs GRE with FatNavs recon on the CUBRIC cluster."
	echo ""
	echo "Usage: cluster_submitGRErecon.sh <input_raw_twix_file>"
	echo ""
	echo "IMPORTANT:"
	echo "Before running you will need to make sure that you set"
	echo "the paths to the retroMoCoBox, the MIRT toolbox and SPM"
	echo "correctly inside this script."
	echo ""
	exit 1;
fi

inputfile=$1

export RETROMOCOBOX_HOME=/home/scedg10/retroMoCoBox
export MIRT_HOME=/home/scedg10/matlab/matlabdownloads/mirt
export SPM_HOME=/cubric/software/spm.versions/spm12
export CLUSTER_LOG_PATH=${PWD}/clusterlogs

# if 'clusterlogs' folder doesn't exist already, create it:
[ -d clusterlogs ] || mkdir clusterlogs

tempfile="$(mktemp)"
echo "#!/bin/bash" > ${tempfile}
echo "#SBATCH --job-name=GREFatNavs" >> ${tempfile}
echo "#SBATCH -p cubric-centos7" >> ${tempfile}
echo "#SBATCH -o ${CLUSTER_LOG_PATH}/GRE_FatNavs_%j.out" >> ${tempfile}
echo "#SBATCH -e ${CLUSTER_LOG_PATH}/GRE_FatNavs_%j.err" >> ${tempfile}
echo "#SBATCH --mem-per-cpu=16000M" >> ${tempfile}
echo "#SBATCH --ntasks 10" >> ${tempfile}
echo "cd ${RETROMOCOBOX_HOME}" >> ${tempfile}
echo "matlab -nodisplay -nodesktop -nosplash -r \"rawDataFile='${PWD}/${inputfile}';swapDims_xyz = [0 0 1];bKeepGRAPPArecon=1;run_SiemensGRErecon_cluster;exit;\"" >> ${tempfile}

sbatch ${tempfile}

rm ${tempfile}
