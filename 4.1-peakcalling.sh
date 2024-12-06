#!/bin/bash

#SBATCH --partition=queue_name       # the requested queue
#SBATCH --nodes=1              # number of nodes to use
#SBATCH --tasks-per-node=1     # for parallel distributed jobs
#SBATCH --cpus-per-task=8      # for multi-threaded jobs
#SBATCH --mem-per-cpu=2G      # in megabytes, unless unit explicitly stated
#SBATCH --time=40:00:00
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

module load py-numpy/1.26.1-4g7o5u4
module load py-macs2/2.2.8-kx3pu5m

export workingdir=your/working/dir

##REMEMBER: set up any directories that the software needs in this script in case
##it is unable to do so itself

mkdir peaks

#################################################################################
# Main CMD
#################################################################################

# Create a statistical model to define what mapped reads form a peak

# Note1: macs2 -g is effective or mappable genome size, it will vary depending on
# species of interest

# Note2: macs2 --keep-dup: default behaviour is to only leave one optical
# duplicate. If picard is run, then --keep-dup should be kept as all

# Note3: macs2 --broad: add flag if analysing histone peaks

# Note4: macs2 -B will create bedGraph files for fragment pileup if building
# lambda model

declare -a files

for file in $workingdir/bowtie/*
do
        # Extrating the sample names whilst also ensuring it is a sorted.bam file

        if [[ $file == *.sorted.bam ]]
        then
                files+=("$(basename ${file::-11})")
        fi

done

# Ignore cases

shopt -s nocasematch

# Defining the genotypes

GENOTYPES=("Col_0" \
        "TCP4")

declare -a SAMPLES

# separating beween input and non-input samples

for file in ${files[@]}
do

        if [[ ${file} == *_ab_* || ${file} == *_noab_* ]]
        then
                SAMPLES+=(${file})

        elif [[ $file == *input* ]]
        then
                INPUTS+=(${file})
        fi
done

echo "SAMPLES" ${SAMPLES[@]}
echo "INPUTS" ${INPUTS[@]}

echo "============================="

for genotype in ${GENOTYPES[@]}
do
        # need to state the number of replicates for the loop to work
        # in this case, there are 2 reps so while loop is set to <3
        # rep is the loop counter

        rep=1
        count=3

        echo ${genotype} "samples"
        echo "============================="


        while [ ${rep} -lt ${count} ]
        do
                for sample in ${SAMPLES[@]}
                do
                        # match the sample to correct genotype and rep no

                        if [[ ${sample} == ${genotype}_*_batch${rep}* ]]
                        then
                                echo ${sample}

                                # specify whether sample is ip or -ve control

                                if [[ ${sample} == *_ab_* ]]
                                then
                                        sampletype=("ab")

                                elif [[ ${sample} == *_noab_* ]]
                                then

                                        sampletype=("noab")
                                fi

                                for i in ${INPUTS[@]}
                                do
                                        if [[ ${i} == ${genotype}*batch${rep}* ]]
                                        then
                                                input=${i}
                                        fi
                                done


                                # Call peaks for sample ips using corresponding
                                # input as control

                                echo ${sample}".sorted.bam"
                                echo "AGAINST"
                                echo ${genotype} "input" ${rep}
                                echo ${input}".sorted.bam"
                                echo "running..."

                                macs2 callpeak \
                                        -t $workingdir/bowtie/${sample}.sorted.bam \
                                        -c $workingdir/bowtie/${input}.sorted.bam \
                                        -f BAMPE  \
                                        -g 119146348 \
                                        --outdir peaks \
                                        -n ${sample} \
                                        -B

                                echo "complete"
                                echo "======================"


                        fi

                done

        rep=$((${rep} + 1))

        done
        rep=1
done

#################################################################################
# End
#################################################################################