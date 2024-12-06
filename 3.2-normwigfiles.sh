#!/bin/bash

#SBATCH --partition=queue_name       # the requested queue
#SBATCH --nodes=1              # number of nodes to use
#SBATCH --tasks-per-node=1     # for parallel distributed jobs
#SBATCH --cpus-per-task=8      # for multi-threaded jobs
#SBATCH --mem-per-cpu=5G      # in megabytes, unless unit explicitly stated
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

module load samtools/1.10
module load bedtools/2.29.1
module load deeptools/3.5.1

## point to the working directory

export workingdir=your/working/dir

mkdir wigs/wigsnorm

##REMEMBER: set up any directories that the software needs in this script in case
##it is unable to do so itself

#################################################################################
# Main CMD
#################################################################################

declare -a files

for file in $workingdir/bowtie/*sorted.bam
do
        # Arbitrarilly taking the R1 lanes to extrate the sample name

        files+=("$(basename ${file::-11})")

done


echo "bamCoverage = running..."

for file in ${files[@]}
do
       echo "============================="
       echo ${file} "= running..."

## Use the aligned and indexed .sorted.bam files to create normalised bw files
## with the Bin Per Million option (similar to TPM in RNA-seq). Also set smooth
## to a size greater than the bin size to get continuous plot and centre reads for
## sharper signal in enriched regions

       bamCoverage \
               -b bowtie/${file}.sorted.bam \
               -o $workingdir/wigs/wigsnorm/${file}.bw \
               --binSize 20 \
               --normalizeUsing BPM \
               --smoothLength 60 \
               --extendReads 150 \
               --centerReads \
               -p 6 2> $workingdir/wigs/wigsnorm/${file}_bamCoverage.log

       echo ${file} "complete"
       echo "============================="
done

echo "bamCoverage === complete"
echo "============================="

## Same again but this time we are also normalising against the input

# Ignore cases

shopt -s nocasematch

# Defining the genotypes

GENOTYPES=("Col_0" \
        "TCP4")

declare -a INPUTS
declare -a SAMPLES

# Sort betweeen samples and inputs

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

echo "Inputs = " ${INPUTS[@]}
echo "Samples =" ${SAMPLES[@]}

echo "bamCompare = running..."
echo "============================="

for genotype in ${GENOTYPES[@]}
do

# Replicate counter

        rep=1
        count=3


        echo ${genotype} "samples"
        echo "============================="


        while [ ${rep} -lt ${count} ]
        do
                for sample in ${SAMPLES[@]}
                do
                        # Match the sample to correct genotype and rep no

                        if [[ ${sample} == ${genotype}_*_batch${rep}* ]]
                        then

                                for i in ${INPUTS[@]}
                                do
                                        if [[ ${i} == ${genotype}*batch${rep}* ]]
                                        then
                                                input=${i}
                                        fi
                                done


                                # Normalise sample using corresponding
                                # input

                                echo ${sample}".sorted.bam"
                                echo "AGAINST"
                                echo ${input}".sorted.bam"
                                echo "running..."

                        bamCompare \
                                -b1 bowtie/${sample}.sorted.bam \
                                -b2 bowtie/${input}.sorted.bam \
                                -o $workingdir/wigs/wigsnorm/${sample}_vs_input.bw \
                                --binSize 20 \
                                --normalizeUsing BPM \
                                --smoothLength 60 \
                                --extendReads 150 \
                                --centerReads \
                                --scaleFactorsMethod None
                                -p 6 2> $workingdir/wigs/wigsnorm/${sample}_vs_input_bamCoverage.log


                        fi

                done

        rep=$((${rep} + 1))

        done
        rep=1

done

echo "bamCompare = complete"
echo "============================="

#################################################################################
# End
#################################################################################