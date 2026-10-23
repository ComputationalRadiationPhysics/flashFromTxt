module purge

module rome

module load gcc/7.3.0
module load openmpi
module load zlib
module load gnuplot/6.0.4
#module load netcdf/4.7.4
alias sq='squeue -u $USER'

#module load hdf5-parallel/1.8.20
export HDF5_PREFIX="$HOME/lib/gcc7.3.0-openmpi4.1.8/hdf5-1.14.6-parallel-mpif90"
export HYPRE_PREFIX="$HOME/lib/gcc7.3.0-openmpi4.1.8/hypre-2.30.0-mpif90"

export PATH="$HDF5_PREFIX/bin:$PATH"
export LD_LIBRARY_PATH="$HDF5_PREFIX/lib:$HYPRE_PREFIX/lib:$LD_LIBRARY_PATH"
export CMAKE_PREFIX_PATH="$HDF5_PREFIX:$HYPRE_PREFIX:$CMAKE_PREFIX_PATH"
export PKG_CONFIG_PATH="$HDF5_PREFIX/lib/pkgconfig:$HYPRE_PREFIX/lib/pkgconfig:$PKG_CONFIG_PATH"
#module load idl/7.1
module load lapack
module load python


# "tbg" default options #######################################################
#   - SLURM (sbatch)
#   - "defq" queue

# allocate an interactive shell for one hour
#   getNode 2  # allocates two interactive nodes (default: 1)
function getNode() {
    if [ -z "$1" ] ; then
        numNodes=1
    else
        numNodes=$1
    fi
    srun --time=1:00:00 --nodes=$numNodes --ntasks-per-node=2 --cpus-per-task=20  --mem=360000 -p cpu-rome --pty bash
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
    srun --time=1:00:00 --ntasks-per-node=$(($numDevices)) --cpus-per-task=$((20 * $numDevices)) --mem=$((180000 * numDevices)) -p cpu-rome --pty bash
}

function get3HoursDevice() {
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
    srun --time=3:00:00 --ntasks-per-node=$(($numDevices)) --cpus-per-task=$((20 * $numDevices)) --mem=$((180000 * numDevices)) -p cpu-rome --pty bash
}



function get24HoursDevice() {
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
    srun --time=24:00:00 --ntasks-per-node=$(($numDevices)) --cpus-per-task=$((20 * $numDevices)) --mem=$((180000 * numDevices)) -p cpu-rome --pty bash
}