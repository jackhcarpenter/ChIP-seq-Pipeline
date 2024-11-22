#!/bin/bash

#SBATCH --partition=queue_name       # the requested queue
#SBATCH --nodes=1              # number of nodes to use
#SBATCH --tasks-per-node=1     # for parallel distributed jobs
#SBATCH --cpus-per-task=6      # for multi-threaded jobs
#SBATCH --mem-per-cpu=2G      # in megabytes, unless unit explicitly stated
#SBATCH --error=logs/%J.err         # redirect stderr to this file
#SBATCH --output=logs/%J.out        # redirect stdout to this file
#SBATCH --mail-user=your.email@host     # email
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

export workingdir=your/working/dir

echo "working dir =" $workingdir

mkdir fastp

export exportdir=your/working/dir/fastp

echo "export dir =" $exportdir

##REMEMBER: set up any directories that the software needs in this script in case
##it is unable to do so itself

#################################################################################
# Main CMDs
#################################################################################

# Loop variables

# Creating an array containing one instance of each sample ID

declare -a files

for file in $workingdir/merged/*
do
        # Arbitrarilly taking the R1 lanes to extrate the sample name whilst
        # also ensuring it is a fastq file
        if [[ $file == *R1.fastq.gz ]]
        then
                files+=("$(basename ${file::-12})")
        fi

done

echo ${files}

# Trim low quality reads, remove adapters, and poly Gs

echo "RUNNING fastp"

for file in ${files[@]}
do
        echo ${file} "= running"

        fastp \
            -i $workingdir/merged/${file}_R1.fastq.gz \
                -I $workingdir/merged/${file}_R2.fastq.gz \
            --detect_adapter_for_pe \
            --trim_poly_g \
            --correction \
            -o $exportdir/${file}_R1.fastp \
            -O $exportdir/${file}_R2.fastp

        echo ${file} "= complete"

done


echo "fastp COMPLETE"
echo "============================="
#################################################################################
# End
#################################################################################