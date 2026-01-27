#!/bin/bash

## Purpose: Turbomole JOB example script for bwHPC, such as bw{For,Uni}Cluster
##          for   S I N G L E   N O D E   runs  on $TMPDIR  O N L Y
##
## History:
## -------- 
##   original: C. Mosch, UUlm, Version 6.2 2012-05-07
##  last mod.: R. Barthel, KIT, 2020-04-14
##        $Id: $

## How to read comments:
## ---------------------
## 1. Every line starting with '##' is an official comment.
## 2. Any command given in the comment section starts with a '$'.
##    To execute any command from the comment section, please remove '$' first
## 3. 

## Notes on modification of this script:
## -------------------------------------
## 1. FOR MASS JOB SUBMITS, PLEASE REMOVE ALL NOT REQUIRED COMMENTS IN THIS FILE, via:
##      $ sed -e '/^#/d' bwHPC_turbomole_single-node_tmpdir_example.sh > my_TM_script.sh 
##
## 2. PLEASE do not generate LOCAL COPIES of this JOB example script. Be up-to-date and 
##    always copy the JOB example script from the Turbomole example directory.

## Quickstart:
## -----------
## 1.  Login to bwHPC Cluster, e.g. bwUniCluster via:
##       $ ssh <username>@bwunicluster.scc.kit.edu
## 2.  Create a workspace dir for temporarily storing the results of your jobs:
##       $ ws_allocate calc_repo 30; cd $(ws_find calc_repo)
##     Please be aware that the data stored there is only of temporary nature.
## 3.  For each job, create a new subdirectory in your workspace and change to it:
##       $ cd $(ws_find calc_repo); mkdir my_first_job; cd my_first_job
## 4.  Copy the Turbomole example input and queuing system script to the new dir:
##       $ module load chem/turbomole
##       $ cp -v $TURBOMOLE_EXA_DIR/* ./
## 5.  Submit the job in the directory 'my_first_job' to the developer queue via command:
##       $ sbatch -p dev_single bwHPC_turbomole_single-node_tmpdir_example.sh
##     The sbatch command returns the job-ID of your new job.
##     NOTE: For your real jobs, use have to specify a dffferent job queue
## 6.  To monitor or cancel your job, please read:
##        https://wiki.bwhpc.de/e/bwUniCluster_2.0_Batch_Queues
## 7.  The stdout and stderr of running and finished jobs can be found in the
##     submit dir (e.g. bwUniCluster -> slurm_123456.out) of the job.
## 8.  A result tgz file (e.g. bwUniCluster -> job_uc2_123456.tgz)
##     will be copied back to the submit DIR at the end of the job.

## Turbomole and FILE SYSTEM:
## --------------------------
## 1.  Please do NOT calculate in your $HOME or submit directory. Use a local
##     temporary work DIR $TMP, $TMPDIR or $SCRATCH for your calculations.
## 2.  Please do NOT copy files around with 'scp', use 'cp' or 'rsync' instead. Your submit
##     directory is mounted on every compute node.
## 3.  Please save disk space: Delete unnecessary files and compress result
##     files _before_ copying them back to your submit directory.

## USE of SBATCH and its environment variables:
## --------------------------------------------
## *   For more information read https://wiki.bwhpc.de/e/BwUniCluster_2.0_Batch_Queues
## Note  :
##   - all '#SBATCH -?' are equivalent to 'sbatch -?'
##   - for more options, see 'man sbatch'
##   - '# SBATCH' is a comment, while '#SBATCH' is a sbatch statement

#########################################
## E D I T  here your SBATCH resources ##
#########################################
#SBATCH --job-name TM_test_mpi
#SBATCH --ntasks=2 --nodes=1
#SBATCH --mem-per-cpu=2000
#SBATCH --time=00:10:00

#########################################
## Start of job shell script           ##
#########################################
echo -e "\n### Setting up shell environment ...\n"
unset LANG; unset LC_CTYPE; export MKL_NUM_THREADS=1; export OMP_NUM_THREADS=1
export QS_USER=${SLURM_JOB_USER:=$(logname)}
export QS_JOBID=${SLURM_JOB_ID:=$(date +%s)}
export QS_SUBMITDIR=${SLURM_SUBMIT_DIR:=$(pwd)}      
export QS_JOBNAME=${SLURM_JOB_NAME:=$(basename "$0")}
export QS_NPROCS=${SLURM_NPROCS:=2}
export QS_NNODES=${SLURM_NNODES:=1}
## Remarks:
## * The LANG + LC_CTYPE commands ensure that floats are written with a dot ".".
## * MKL_NUM_THREADS=1 and OMP_NUM_THREADS=1 disable all automatic threading
##   mechanisms in MKL (Intel Math Kernel Library). So the number of cores are
##   under sole control of Turbomole (see PARNODES below).
## * Default values are assigned if queueing system variables
##   are not defined (e.g. if accidentally running the script without sbatch).

echo -e "\n### Printing basic job infos to stdout ...\n"
echo "START_TIME           = $(date +'%Y-%m-%d %H:%M:%S')"
echo "HOSTNAME             = ${HOSTNAME}"
echo "JOB USER             = ${QS_USER}"
echo "JOB NAME             = ${QS_JOBNAME}"
echo "JOB ID               = ${QS_JOBID}"
echo "Submit directory     = ${QS_SUBMITDIR}"
echo "Number of CPUs       = ${QS_NPROCS}"
echo "Number of nodes      = ${QS_NNODES}"

echo -e "\n### Creating RUN directory on local disc and changing to it ...\n"
## Since we support parallelism on one node only (up to 16 cores),
## we always use $TMPDIR as run directory.
## 
## NEVER EVER calculate in your home directory.
##

## Remarks about env variable $TMPDIR
## * bwUniCluster:           TMPDIR=/scratch/slurm_tmpdir/Job_<Job ID>
## * bwForCluster Chemistry: a) diskless nodes: TMPDIR=/tmp/<username>_job_<JOB ID>
##                           b) nodes w disks : TMPDIR=/scratch/<username>_job_<JOB ID>

## Check $runDIR and generate if neccessary
runDIR=${TMPDIR}
##
## a) Avoid that $runDIR points to root dir of /scratch or /tmp
if [[ "${runDIR}" == "/scratch" ]] ; then runDIR="/scratch/$USER/Job_${QS_JOBID}" ; fi
if [[ "${runDIR}" == "/tmp"     ]] ; then runDIR="/tmp/$USER/Job_${QS_JOBID}" ; fi
##
## b) Setup variables for tar archive
basename_runDIR="$(basename $runDIR)"
dirname_runDIR="$(dirname $runDIR)"
tararchive=$(echo "${basename_runDIR}" | sed -e "s/$USER//;s/[j,J]ob/job_${CLUSTER}/;s/^_//")
echo "runDIR               = ${runDIR}"
echo "Basename of runDIR   = ${basename_runDIR}"
echo "Dirname  of runDIR   = ${dirname_runDIR}"
echo "Tar archive name     = ${tararchive}"
##
## c) Generate $runDIR
mkdir -vp ${runDIR}
##
## d) Change to $runDIR 
cd ${runDIR}


echo -e "\n### Loading modules and defining environment:\n"

## Unload first all loaded turbomole versions
module unload chem/turbomole


## Define TURBOMOLE_MODE+PARA_ARCH (must be done before 'module load' command):
export TURBOMOLE_MODE="compute"


if [ "$QS_NPROCS" -gt 1 ] ; then
##   Please be aware that some of the Tubomole commands are shared memory
##   parallelized (PARA_ARCH=SMP) and some of the commands are parallelized
##   via message passing intereface (PARA_ARCH=MPI). Please read
##   chapter "Parallel runs" in "Tutorial_7-3-1.pdf" and
##   chapter "Parallel Runs" in "Documentation.pdf" for
##   documentation about which Turbomole command offers which parallelization.
##   Usually one selects "SMP" for jobs, that run entirely within one node
##   and "MPI" for jobs, that span multiple nodes. Nevertheless in some
##   cases (depending on job input) the MPI version can be faster than
##   the SMP version, even on one node. Therefore before submitting mass jobs
##   you should always benchmark the speed of your jobs by modifying PARA_ARCH
##   and PARNODES for the same job and comparing the reached speed up.
  
  export PARA_ARCH="SMP" 

##   Please make sure that the number of cores requested via queueing system
##   matches the number of cores stored in PARNODES.
  export PARNODES="$QS_NPROCS"

else

  unset PARNODES
  unset PARA_ARCH

fi

module load chem/turbomole  ## this loads default Turbomole version

if [ -z "${TURBOMOLE_VERSION}" ] ; then
  echo "ERROR: Could not load module 'chem/turbomole'."
  exit -1
fi
echo -e "Loaded module        = ${LOADEDMODULES//:/ }\n"

## Check $TURBOTMPDIR
## a) for PARA_ARCH=SMP or PARA_ARCH="" -> $TURBOTMPDIR must be equal $runDIR
## b) for PARA_ARCH=MPI -> $TURBOTMPDIR must be equal $runDIR + token
## c) for PARA_ARCH=MPI and multinode jobs -> $TURBOTMPDIR must be emptied
if [ "${PARA_ARCH}" == "MPI" ] ; then
   if [ ${QS_} > 1 ] ; then
      ## Multinode MPI jobs can only be run on a global filesystem and there $TURBOTMPDIR should be deactivated
      export TURBOTMPDIR="" 
   else
      ## TURBOTMPDIR will set $tmpdir in control, which for MPI runs is the token that generates the subdir directory names
      ## e.g. $TURBOTMPDIR=/scratch/run-dir -> /scratch/run-dir-001/dens
      export TURBOTMPDIR="${runDIR}/run-dir"
   fi
else
   export TURBOTMPDIR="${runDIR}"
fi
## 
echo "TURBOMOLE_VERSION    = $TURBOMOLE_VERSION"
echo "TURBOTMPDIR          = $TURBOTMPDIR (\$runDIR=$runDIR)"
echo "TURBOMOLE_MODE       = $TURBOMOLE_MODE"
echo "PARA_ARCH            = $PARA_ARCH"
echo "PARNODES             = $PARNODES"


echo -e "\n### Copying input files for job (if required):\n"

## In this section the Turbomole input files '{coord,*basis,control,mos}'
## are copied from your submit directory to your local temporary working dir.
##
## The input files can be created via command line tool 'define' (as described
## in 'Tutorial_7-4.pdf') or via graphical user interface 'TmoleX' (as described
## in 'Tutorial-tmolex-4-5.pdf'). 
## If you are using TmoleX, please only construct and save jobs ('Save'
## button) via 'TmoleX', but do NEVER use the 'Run'-job buttons (since they do not
## work properly). Furthermore you can extract the Turbomole commands (e.g.
## 'dscf', 'jobex', ...) of your job from file 'start-job' saved by TmoleX
## and paste these commands below into this queueing system script.

cp -v $QS_SUBMITDIR/{coord,*basis,control,mos}         $runDIR/
# cp -v $TURBOMOLE_EXA_DIR/{coord,*basis,control,mos}      $runDIR/
# cp -v $HOME/some_fixed_path/{coord,*basis,control,mos}   $runDIR/

## TURBOMOLE_EXA_DIR is a dir with some simple Turbomole example input files.
## QS_SUBMITDIR is the dir in which you did type 'sbatch ...'.


echo -e "\n### Performing Turbomole calculations ...\n"
## Here you usually add a series of Turbomole commands calling
## different Turbomole sub-programs, e.g. 'dscf', 'escf', 'grad', 'ricc2', ...
## In case of parallel jobs, please check that the value of PARA_ARCH (see
## above) is appropriate for the called sub-programs.
##
## If you use MS-Windows please make sure to convert <CR><LF> to <CR> in all
## files before submitting the job (see Linux command 'dos2unix').

time dscf > dscf.out 2>&1



echo -e "\n### Cleaning up files ... removing unnecessary scratch files ...\n"
rm -vf dens errvec fock oldfock slave*
sleep 10 
## Remarks:
## * PLEASE remove all not required files in this script.
## * The "sleep 10" command helps to prevent stale (dead) nfs files.



echo -e "\n### Compressing results and copying back result archive ...\n"
cd "${dirname_runDIR}"
mkdir -vp "${QS_SUBMITDIR}" ## if user has deleted or moved the submit dir
echo "Creating tgz-file '${QS_SUBMITDIR}/${tararchive}.tgz' ..."
tar -zcvf "${QS_SUBMITDIR}/${tararchive}.tgz" "${basename_runDIR}"
## Remarks:
## * The resulting tgz file is copied back to the submit directory.
##   The name of the tgz file looks similar too "job_uc1_123456.tgz"



echo -e "\n### Cleanup of $runDIR is done by SYSTEM automatically\n"
echo "END_TIME             = $(date +'%Y-%m-%d %H:%M:%S')"

