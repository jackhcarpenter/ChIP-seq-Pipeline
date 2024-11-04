#!/bin/bash

#SBATCH --partition=QUEUE_NAME       # the requested queue
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

module load fastqc/v0.11.9
module load fastp/v0.20     ## NOTE: if you are unable to disable paralisation
module load multiqc/v/1.9   ## then you must run fastp and multiqc seperatly as
                            ## they use conflicting versions of python

export workingdir=/your/working/dir

##REMEMBER: set up any directories that the software needs in this script in case 
##it is unable to do so itself

#################################################################################
# Main CMD
#################################################################################

# Creating an array containing one instance of each sample ID

# need to be able to separate by lanes
lanes=("L001" \
        "L002")

declare -a files

for file in $workingdir/fastp/*
do
        # Arbitrarilly taking the R1 lanes to extrate the sample name whilst
        # also ensuring it is a fastq file
        if [[ $file == *R1*.fastq.gz ]]
        then
                files+=("$(basename ${file::-14})")
        fi

done

## perform fastqc on the trimmed PE data

for lane in ${lanes[@]}
do

        for file in ${files[@]}
        do
                echo ${file} "running"

                fastqc $workingdir/${file}_${lane}_R1.fastp
                fastqc $workingdir/${file}_${lane}_R2.fastp

                echo ${file} "complete"

        done

        ## summarise the QC data of all reads
        multiqc -i "TCP4_ChIP_LANE_"${lane} $workingdir/

        echo "multiqc for "${lane}" complete"

done

echo ${"QC complete"}
echo ${"============================="}

#################################################################################
# End
#################################################################################