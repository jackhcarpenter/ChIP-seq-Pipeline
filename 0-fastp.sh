#!/bin/bash

#SBATCH --partition=jumbo       # the requested queue
#SBATCH --nodes=1              # number of nodes to use
#SBATCH --tasks-per-node=1     # for parallel distributed jobs
#SBATCH --cpus-per-task=4      # for multi-threaded jobs
#SBATCH --mem-per-cpu=16G      # in megabytes, unless unit explicitly stated
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

echo "working dir =" $workingdir

mkdir fastp

export exportdir=/mnt/scratch/c1831460/ChIP/fastp

echo "export dir =" $exportdir

##REMEMBER: set up any directories that the software needs in this script in case
##it is unable to do so itself

#################################################################################
# Main CMDs
#################################################################################

# Loop variables

# Creating an array containing one instance of each sample ID

#need to be able to separate by lanes
lanes=("L001" \
        "L002")

declare -a files

for file in $workingdir/RAW_DATA/S349_NovaSeq_BHLVN2DRX5/fastq/*
do
        # Arbitrarilly taking the R1 lanes to extrate the sample name whilst
        # also ensuring it is a fastq file
        if [[ $file == *R1*.fastq.gz ]]
        then
                files+=("$(basename ${file::-21})")
        fi

done

echo ${files}

# Trim low quality reads, remove adapters, and poly Gs

echo "RUNNING fastp"

for lane in ${lanes[@]}
do

        for file in ${files[@]}
        do
                echo ${file} "= running"

                fastp \
                    -i $workingdir/RAW_DATA/S349_NovaSeq_BHLVN2DRX5/fastq/${file}_${lane}_R1_001.fastq.gz \
                        -I $workingdir/RAW_DATA/S349_NovaSeq_BHLVN2DRX5/fastq/${file}_${lane}_R2_001.fastq.gz \
                    --detect_adapter_for_pe \
                    --trim_poly_g \
                    --correction \
                    -o $exportdir/${file}_${lane}_R1.fastp \
                    -O $exportdir/${file}_${lane}_R2.fastp

                echo ${file} "= complete"

        done

        echo ${lane} "complete"
done

echo "fastp COMPLETE"
echo "============================="
#################################################################################
# End
#################################################################################