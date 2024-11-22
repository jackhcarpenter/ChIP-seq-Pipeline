#!/bin/bash

#SBATCH --partition=queue_name       # the requested queue
#SBATCH --nodes=1              # number of nodes to use
#SBATCH --tasks-per-node=1     # for parallel distributed jobs
#SBATCH --cpus-per-task=4      # for multi-threaded jobs
#SBATCH --mem-per-cpu=2G      # in megabytes, unless unit explicitly stated
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

module load fastqc/
module load multiqc/

export workingdir=your/working/dir

##REMEMBER: set up any directories that the software needs in this script in case
##it is unable to do so itself

mkdir fastqc

#################################################################################
# Main CMD
#################################################################################

# Creating an array containing one instance of each sample ID

declare -a files

for file in $workingdir/fastp/*
do
        # Arbitrarilly taking the R1 lanes to extrate the sample name whilst
        # also ensuring it is a fastq file
        if [[ $file == *R1.fastp ]]
        then
                files+=("$(basename ${file::-9})")
        fi

done

## perform fastqc on the trimmed PE data

    #    for file in ${files[@]}
   #     do
  #              echo ${file} "running"

 #               fastqc $workingdir/fastp/${file}_R1.fastp \
#                       -o $workingdir/fastqc
 #               fastqc $workingdir/fastp/${file}_R2.fastp \
#                       -o $workingdir/fastqc

#                echo ${file} "complete"

#        done

        ## summarise the QC data of all reads

## summarise the QC data of all reads

#mv fastp/*_fastqc.* fastqc/

multiqc -i "TCP4_ChIP" fastqc/ \
        --ignore unmerged \
        -o $workingdir/fastqc

echo "MultiQC complete"
echo "============================="

#################################################################################
# End
#################################################################################
