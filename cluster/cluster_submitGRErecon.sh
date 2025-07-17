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
	echo "multi-echo data. Please make sure you have enough quota to handle this!"
	echo ""
	echo "IMPORTANT:"
	echo "Before running you will need to make sure that you set"
	echo "the RETROMOCOBOX_HOME environment variable (and also ASPIRE_HOME"
	echo "if you are using ASPIRE)."
	echo "Alternatively, these can be set by manually editing this script and"
	echo "specifying them below this text in the .sh file itself."
	exit 1;
fi

export RETROMOCOBOX_HOME=/home/scedg10/retroMoCoBox
export ASPIRE_HOME=/home/scedg10/matlab_cubric/ASPIRE

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

if [[ -z "${ASPIRE_HOME}" ]]; then
    echo ""
    echo "ERROR: before running this script please ensure your have set"
    echo "the ASPIRE_HOME variable to point to your local version"
    echo "of the ASPIRE matlab code. This can be set to a dummy value"
    echo "if your data is not multi-echo."
    echo ""
    exit 1;
fi


export SPM_HOME=/cubric/software/spm.versions/spm12

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
echo "matlab -nodisplay -nodesktop -nosplash -r \"rawDataFile='${inputfile}';outRoot='${outRoot}';swapDims_xyz = [0 0 1];bKeepGRAPPArecon=0;parpoolSize = 6;script_pwd = '${PWD}';ASPIRE_HOME='${ASPIRE_HOME}';run_SiemensGRErecon_cluster;exit;\"" >> ${tempfile}
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
