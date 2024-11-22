#!/bin/bash

#SBATCH --partition=queue_name       # the requested queue
#SBATCH --nodes=1              # number of nodes to use
#SBATCH --tasks-per-node=1     # for parallel distributed jobs
#SBATCH --cpus-per-task=4      # for multi-threaded jobs
#SBATCH --mem-per-cpu=4G      # in megabytes, unless unit explicitly stated
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

module load bowtie2/v2.4.1
module load samtools/1.10
module load bamtools/v2.5.1

## point to the directory containing the reference genome where sequences will be
## mapped
export refdir=your/working/dir/At_reference_genome

## point to the working directory
export workingdir=your/working/dir

##REMEMBER: set up any directories that the software needs in this script in case
##it is unable to do so itself

mkdir bowtie/

#################################################################################
# Main CMD
#################################################################################

## List of sequences to map to indexed reference genome

declare -a files

for file in $workingdir/fastp/*
do
        # Arbitrarilly taking the R1 lanes to extrate the sample name

        if [[ $file == *R1.fastp ]]
        then
                files+=("$(basename ${file::-9})")

        else
                echo $file "is in wrong format"
        fi

done


for file in ${files[@]}
do

        echo "============================="
        echo ${file} "= running mapping"

        ## Map forward and reverse reads to the indexed referenced genome

        bowtie2 \
        --maxins 500 \
        --fr \
        -p 4 \
        -x $refdir/Arabidopsis_thaliana.TAIR10.59.gtf.gz \
        -1 $workingdir/fastp/${file}_R1.fastp \
        -2 $workingdir/fastp/${file}_R2.fastp \
        -S $workingdir/bowtie/${file}.sam

        ## Compress the aligned sam files to bam files
        samtools view \
        -b $workingdir/bowtie/${file}.sam \
        > $workingdir/bowtie/${file}.bam \

        echo ${file} "= mapping complete"

        ## Organise mapped reads and index them for faster access during
        ## downstream processing

        echo ${file} "= sorting"

        samtools sort \
        -@ ${SLURM_CPUS_PER_TASK} \
        -o $workingdir/bowtie/${file}.sorted.bam \
        $workingdir/bowtie/${file}.bam

        samtools index \
        $workingdir/bowtie/${file}.sorted.bam

        echo ${file} "= sorting complete"

        ## Run some stats on aligned/mapped reads
        bamtools stats \
        -in $workingdir/bowtie/${file}.sorted.bam \
        > $workingdir/bowtie/${file}.sorted.stats.txt
        ## look at the stats files to quality check the data

        echo ${file} "= stats complete"
        echo "============================="

done
