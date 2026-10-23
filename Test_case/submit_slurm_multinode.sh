#!/bin/bash -l
#SBATCH --job-name=4P16_parallel_2cpus
#SBATCH --partition=cpu-rome
#SBATCH --time=96:00:00
#SBATCH --nodes=1                 # number of nodes
#SBATCH --ntasks-per-node=2     # cpu-rome: 128 cores / 512 GB per node
#SBATCH --cpus-per-task=1         # pure MPI: one core per rank
#SBATCH --mem-per-cpu=32G          # 128 x 3G = 384G < 512G node; ~1 GB/rank needed with maxblocks=1000
#SBATCH --output=flash_%j.out
#SBATCH --mail-type=ALL
##SBATCH --mail-user=<your.name@hzdr.de>

# --- same environment as the build, on every node -------------------------
source ~/FLASH_rome.sh
export OMP_NUM_THREADS=1               # no hidden threading from Hypre/BLAS
export HDF5_USE_FILE_LOCKING=FALSE     # parallel HDF5 1.14 on shared FS (GPFS/NFS)

cd "$SLURM_SUBMIT_DIR"

echo "Job $SLURM_JOB_ID on $SLURM_JOB_NUM_NODES nodes: $SLURM_JOB_NODELIST"
echo "Total MPI ranks: $SLURM_NTASKS"
echo "--- libraries flash4 will use (must all be the same OpenMPI) ---"
ldd ./flash4 | grep -Ei 'libmpi|hdf5|HYPRE'
echo "------------------------------------------------------------------"

# Let OpenMPI take the rank count and host list from Slurm.
# -x forwards the environment to the ranks on the remote nodes.
mpirun -np "$SLURM_NTASKS" \
       --map-by core --bind-to core \
       -x PATH -x LD_LIBRARY_PATH -x OMP_NUM_THREADS -x HDF5_USE_FILE_LOCKING \
       ./flash4

# Alternative if your OpenMPI is built with Slurm PMIx support (check: srun --mpi=list):
# srun --mpi=pmix ./flash4

echo 'done'