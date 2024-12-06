#!/bin/bash

#SBATCH --partition=queue_name      # the requested queue
#SBATCH --nodes=1              # number of nodes to use
#SBATCH --tasks-per-node=1     # for parallel distributed jobs
#SBATCH --cpus-per-task=6      # for multi-threaded jobs
#SBATCH --mem-per-cpu=1G      # in megabytes, unless unit explicitly stated
#SBATCH --time=20:00:00
#SBATCH --error=logs/%J.err         # redirect stderr to this file
#SBATCH --output=logs/%J.out        # redirect stdout to this file
#SBATCH --mail-user=your.email@host # email address used for event notification
#SBATCH --mail-type=BEGIN,END,FAIL # email on job start, end, and/or failure

#################################################################################
# Print Slurm Parameters to Console
#################################################################################

echo "Usable Environment Variables:"
echo "============================="
echo "hostname=$(hostname)"
echo \$SLURM_JOB_ID=${SLURM_JOB_ID}
echo \$SLURM_NTASKS=${SLURM_NTASKS}
echo \$SLURM_NTASKS_PER_NODE=${SLURM_NTASKS_PER_NODE}
echo \$SLURM_CPUS_PER_TASK=${SLURM_CPUS_PER_TASK}
echo \$SLURM_JOB_CPUS_PER_NODE=${SLURM_JOB_CPUS_PER_NODE}
echo \$SLURM_MEM_PER_CPU=${SLURM_MEM_PER_CPU}

#################################################################################
# Modulels to Load and Setup
#################################################################################

module load gcc/9.5.0
module load deeptools/3.5.1

## Point to reference genome

export refdir=your/working/dir/At_reference_genome

## point to the working directory

export workingdir=your/working/dir/ChIP

mkdir plots/

##REMEMBER: set up any directories that the software needs in this script in case
##it is unable to do so itself

#################################################################################
# Main CMD
#################################################################################

echo "============================="
echo "computeMatrix running..."

#ls $workingdir/wigs/wigsnorm/TCP4_AB_Batch?_S?.bw

# Taking normalised bw files and creating a matrix that represents average signal 
# per bin relative to the TSS for all peaks for downstream visualisation

computeMatrix reference-point \
	--referencePoint TSS \
	-b 1000 -a 1000 \
	-R $refdir/Arabidopsis_thaliana.TAIR10.59.gtf \
	-S $workingdir/plots/TCP4_AB_Batch*.bw \
	--skipZeros \
	-o $workingdir/plots/TCP4_matrix_TSS.gz \
	-p 6 \
	--outFileSortedRegions $workingdir/plots/regions_TSS.bed

echo "computeMatrix complete"
echo "============================="

echo "plotProfile running..."

# Making a peak distribution plot

plotProfile \
        -m $workingdir/plots/TCP4_matrix_TSS.gz \
        -out $workingdir/plots/TSS_TCP4_profile.png \
        --perGroup \
        --colors green purple \
        --plotTitle "" --samplesLabel "Rep1" "Rep2" \
        --refPointLabel "TSS" \
        -T "TCP4 read density" \
        -z ""

echo "plotProfile complete"
echo "============================="

echo "plotHeatmap running..."

# Similar distribution plot but as a heatmap

plotHeatmap \
        -m  $workingdir/plots/TCP4_matrix_TSS.gz \
        -out $workingdir/plots/TSS_TCP4_heatmap.png \
        --colorMap RdBu \
        --whatToShow 'heatmap and colorbar' \
        --zMin -4 --zMax 4

echo "plotHeatmap complete..."
echo "============================="

#################################################################################
# End
#################################################################################