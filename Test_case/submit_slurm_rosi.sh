#!/bin/bash
#SBATCH --job-name=4P16_2D_Expansion    # Job name
##SBATCH -A fwk    # Job name
#SBATCH --time=96:00:00               # Time limit hrs:min:sec
##SBATCH -p intel,defq -n 5
##SBATCH -p defq -n 5
#SBATCH -p cpu-rome -n 1         # serial: one MPI rank (mpirun -np 1)
#SBATCH -N 1                    # Run on a single CPU
#SBATCH --mail-type=ALL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mem=64G           # maxblocks=22000 x ~1 MB/block ~ 22 GB + solvers, with margin
##SBATCH --mail-user=<your.name@hzdr.de>     # Where to send mail	
 
echo 'start...'
# same environment as the build
source ~/FLASH_rome.sh
#srun --mpi=pmi2 flash4
#srun --mpi=pmi2 ../code_slurm/exe/ipicls2d < ./input.in 
#mpirun -np 1 flash4
mpirun -np 1 flash4
#srun flash4
echo 'done'
 
 