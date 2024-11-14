#!/bin/bash

#SBATCH --partition=queue_name       # the requested queue
#SBATCH --nodes=1              # number of nodes to use
#SBATCH --tasks-per-node=1     # for parallel distributed jobs
#SBATCH --cpus-per-task=4      # for multi-threaded jobs
#SBATCH --mem-per-cpu=4G      # in megabytes, unless unit explicitly stated
#SBATCH --error=logs/%J.err         # redirect stderr to this file
#SBATCH --output=logs/%J.out        # redirect stdout to this file
#SBATCH --mail-user=your.email@host      # email
#SBATCH --mail-type=BEGIN,END,FAIL      # email on job start, end, and/or failure

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

export workingdir=your/working/dir

mkdir merged

##REMEMBER: set up any directories that the software needs in this script in case
##it is unable to do so itself

#################################################################################
# Main CMDs
#################################################################################

## Creating file name list

list=("Col_0_AB_Batch1_S2" \
        "Col_0_AB_Batch2_S8" \
        "Col_0_INPUT_Batch2_S7" \
        "Col_0_Input_Batch1_S3" \
        "Col_0_NOAB_Batch2_S9" \
        "Col_0_NoAB_Batch1_S1" \
        "TCP4_AB_Batch1_S5" \
        "TCP4_AB_Batch2_S11" \
        "TCP4_INPUT_Batch2_S10" \
        "TCP4_Input_Batch1_S6" \
        "TCP4_NOAB_Batch2_S12" \
        "TCP4_NoAB_Batch1_S4" \
        "Undetermined_S0")

## Merging command

for file in ${list[@]}
do
        # combining lane 1 and lane to for forward reads
        echo "==========================================="
        echo "merging" ${file}"_L001_R1 and "${file}"_L002_R1"
        cat $workingdir/RAW_DATA/S349_NovaSeq_BHLVN2DRX5/fastq/${file}_L001_R1_001.fastq.gz \
                $workingdir/RAW_DATA/S349_NovaSeq_BHLVN2DRX5/fastq/${file}_L002_R1_001.fastq.gz \
                >> $workingdir/merged/${file}_R1.fastq.gz

        # combining lane 1 and lane to for reverse reads
        echo "merging" ${file}"_L001_R2 and "${file}"_L002_R2"
        cat $workingdir/RAW_DATA/S349_NovaSeq_BHLVN2DRX5/fastq/${file}_L001_R2_001.fastq.gz \
                $workingdir/RAW_DATA/S349_NovaSeq_BHLVN2DRX5/fastq/${file}_L002_R2_001.fastq.gz \
                >> $workingdir/merged/${file}_R2.fastq.gz

done
echo "==========================================="
echo "Complete"