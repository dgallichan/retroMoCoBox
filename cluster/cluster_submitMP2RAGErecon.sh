#!/bin/bash

if [ $# -lt 1 ];then
	echo ""
	echo "Runs MP2RAGE with FatNavs recon on the CUBRIC cluster."
	echo ""
	echo "Usage: cluster_submitMP2RAGErecon.sh <full_path_to_input_raw_twix_file> <full_path_to_output_root_folder>"
	echo ""
	echo "If you pass just one input file you will get the outputs in the current working folder."
	echo ""
	echo ""
	echo "IMPORTANT:"
	echo "Before running you will need to make sure that you set"
	echo "the RETROMOCOBOX_HOME environment variable"
	echo ""
	exit 1;
fi

inputfile=$1

if [ $# -lt 2 ];then
	outRoot=${PWD}
    export CLUSTER_LOG_PATH=${PWD}/clusterlogs
else
	outRoot=$2
    export CLUSTER_LOG_PATH=${outRoot}/clusterlogs
fi

if [[ -z "${RETROMOCOBOX_HOME}" ]]; then
    echo ""
    echo "ERROR: before running this script please ensure your have set"
    echo "the RETROMOCOBOX_HOME variable to point to your local version"
    echo "of the retroMoCoBox MATLAB tools."
	echo ""
    exit 1;
fi

#export RETROMOCOBOX_HOME=/home/scedg10/retroMoCoBox
export SPM_HOME=/cubric/software/spm.versions/spm12

# if 'cluster log path' folder doesn't exist already, create it:
[ -d ${CLUSTER_LOG_PATH} ] || mkdir ${CLUSTER_LOG_PATH}

tempfile="$(mktemp)"

echo "#!/bin/bash" > ${tempfile}
echo "#SBATCH --job-name=MP2RAGEFatNavs" >> ${tempfile}
echo "#SBATCH -p cubric-centos7" >> ${tempfile}
echo "#SBATCH -o ${CLUSTER_LOG_PATH}/MP2RAGE_FatNavs_%j.out" >> ${tempfile}
echo "#SBATCH -e ${CLUSTER_LOG_PATH}/MP2RAGE_FatNavs_%j.err" >> ${tempfile}
echo "#SBATCH --mem=160G" >> ${tempfile}
echo "#SBATCH --cpus-per-task 12" >> ${tempfile}
echo "#SBATCH --ntasks 1" >> ${tempfile}
echo "cd ${RETROMOCOBOX_HOME}" >> ${tempfile}
echo "matlab2021a -nodisplay -nodesktop -nosplash -r \"rawDataFile='${inputfile}';outRoot='${outRoot}';swapDims_xyz = [0 0 1]; FatNavRes_mm = 2; script_pwd = '${PWD}';run_SiemensMP2RAGErecon_cluster;exit;\"" >> ${tempfile}

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

