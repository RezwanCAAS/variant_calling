# variant_calling
codes are in GATK_2.sh file
##########################################
#                                        #
#    fastp for cleaning the raw reads    #
#                                        #
##########################################

#!/bin/bash
#SBATCH --array=1-64
#SBATCH --ntasks=6
#SBATCH --cpus-per-task=8
#SBATCH --time=30:00:00
#SBATCH --mem-per-cpu=30G
#SBATCH --job-name=DraF1_64_filtering
#SBATCH --output=DraF1_64_filtering.%j.out
#SBATCH --mail-user=tariqr@kaust.edu.sa

module load fastp/0.23.2

SEQDIR='/ibex/scratch/projects/c2141/dragon-fruit/test/'
VAR=($(ls ${SEQDIR}))
OUTDIR='/ibex/scratch/projects/c2141/dragon-fruit/mapped_reads'
