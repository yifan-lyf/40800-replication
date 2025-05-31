#!/bin/bash
#SBATCH --job-name=jags_stage2_0.5
#SBATCH --output=stage2_0.5_out.txt
#SBATCH --error=stage2_0.5_err.txt
#SBATCH --time=20:00:00
#SBATCH --partition=caslake 
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=48G
#SBATCH --account=macs30123

module load R/4.2.2
module load jags/4.3.1

Rscript stage2_0.5.R