#!/bin/bash

#SBATCH --partition=defq       # the requested queue
#SBATCH --nodes=1              # number of nodes to use
#SBATCH --tasks-per-node=1     # for parallel distributed jobs
#SBATCH --cpus-per-task=4      # for multi-threaded jobs
#SBATCH --mem-per-cpu=4G      # in megabytes, unless unit explicitly stated
#SBATCH --error=logs/%J.err         # redirect stderr to this file
#SBATCH --output=logs/%J.out        # redirect stdout to this file
#SBATCH --mail-user=carpenterj3@cardiff.ac.uk      # email
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

module load fastp/v0.20

export workingdir=/mnt/scratch/c1831460/ChIP

##REMEMBER: set up any directories that the software needs in this script in case
##it is unable to do so itself

#################################################################################
# Main CMD
#################################################################################

## list of samples
list=("Col-O-AB_S2" "Col-O-Input_S3" "Col-O-NoAB_S1" "TCPA-AB_S5" "TCPA-Input_S6"\
        "TCPA-NoAB_S4" "Undetermined_S0")

## perform fastp to remove low quality reads and adaptors 
for i in ${list[@]}

do
        echo ${i}

        fastp \
        -i $workingdir/RAW_DATA/230628_fastqs/${i}_merge_R1.fastq.gz \
        -I $workingdir/RAW_DATA/230628_fastqs/${i}_merge_R2.fastq.gz \
        -o $workingdir/trimmed_reads/${i}_fp1.fastq.gz \
        -O $workingdir/trimmed_reads/${i}_fp2.fastq.gz \
        --adapter_sequence=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
	--adapter_sequence_r2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
	--trim_poly_g \
	--correction

done
