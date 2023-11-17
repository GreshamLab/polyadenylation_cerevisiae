#!/bin/sh
#
#SBATCH	--verbose
#SBATCH	--job-name=snakemake
#SBATCH --mail-type=END
#SBATCH	--mail-user=sz4633@nyu.edu
#SBATCH	--output=snake.%j.out
#SBATCH	--error=snake.%j.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --mem=80000
#SBATCH --time 47:59:59

module purge
module load snakemake/6.12.3
module load star/intel/2.7.11a
module load python/intel/3.8.6

export PYTHONUNBUFFERED=TRUE

echo "SLURM Environment: ${SLURM_JOB_NUM_NODES} Nodes ${SLURM_NTASKS} Tasks ${SLURM_MEM_PER_NODE} Memory/Node"
echo ""

time "$@"