module purge
module load gcc
module load openmpi/2.1.2
module load zlib
module load netcdf/4.7.4

#module load hdf5-parallel/1.8.20
export HDF5_ROOT=$HOME/lib/hdf5
export LD_LIBRARY_PATH=$HDF5_ROOT/lib:$HDF5_ROOT/include:$LD_LIBRARY_PATH

module load hypre/2.10.0b
module load idl/7.1
module load lapack
module load python/3.6.5


# "tbg" default options #######################################################
#   - SLURM (sbatch)
#   - "defq" queue
export TBG_SUBMIT="sbatch"
export TBG_TPLFILE="etc/picongpu/hemera-hzdr/defq.tpl"
# use defq for regular queue and defq_low for low priority queue
export TBG_partition="defq"

# allocate an interactive shell for one hour
#   getNode 2  # allocates two interactive nodes (default: 1)
function getNode() {
    if [ -z "$1" ] ; then
        numNodes=1
    else
        numNodes=$1
    fi
    srun --time=1:00:00 --nodes=$numNodes --ntasks-per-node=2 --cpus-per-task=20  --mem=360000 -p defq --pty bash
}

# allocate an interactive shell for one hour
#   getDevice 2  # allocates two interactive devices (default: 1)
function getDevice() {
    if [ -z "$1" ] ; then
        numDevices=1
    else
        if [ "$1" -gt 2 ] ; then
            echo "The maximal number of devices per node is 2." 1>&2
            return 1
        else
            numDevices=$1
        fi
    fi
    srun --time=1:00:00 --ntasks-per-node=$(($numDevices)) --cpus-per-task=$((20 * $numDevices)) --mem=$((180000 * numDevices)) -p defq --pty bash
}

export SCRATCH=/bigdata/hplsim/external/hirsch95