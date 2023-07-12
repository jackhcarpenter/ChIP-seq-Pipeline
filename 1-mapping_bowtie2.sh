#!/bin/bash

#SBATCH --partition=QUEUE_NAME       # the requested queue
#SBATCH --nodes=1              # number of nodes to use
#SBATCH --tasks-per-node=1     # for parallel distributed jobs
#SBATCH --cpus-per-task=8      # for multi-threaded jobs
#SBATCH --mem-per-cpu=4G      # in megabytes, unless unit explicitly stated
#SBATCH --error=logs/%J.err         # redirect stderr to this file
#SBATCH --output=logs/%J.out        # redirect stdout to this file
#SBATCH --mail-user=USERNAME@INSTITUTIONAL_ADDRESS      # email
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
export refdir=/your/reference/directory

## point to the working directory
export workingdir=/your/working/dir

##REMEMBER: set up any directories that the software needs in this script in case 
##it is unable to do so itself

#################################################################################
# Main CMD
#################################################################################

## Index the genome for quicker access by bowtie2 during alignment
bowtie2-build \
    $refdir/SPECIES_GENOME.dna.toplevel.fa\
    $refdir/SPECIES_GENOME.INDEX.fa

## List of sequences to map to indexed reference genome
list=(
        "sample_input1" "sample_inputn" "sample_ip1" "sample_ipn" "sample_neg1"\
        "sample_negn" "control_input1" "control_inputn" "control_ip1" "control_ipn"\
        "control_neg1" "control_negn")

## Map forward and reverse reads to the indexed referenced genome
for i in ${list[@]}
do
        
        echo ${i}

        ## Align the sequences to the genome
        bowtie2 \
        --maxins 500 \
        --fr \
        -p 8 \
        -x $refdir/SPECIES.GENOME.INDEX.fa \
        -1 $workingdir/trimmed_reads${i}_fp1.fastq.gz \
        -2 $workingdir/trimmed_reads/${i}_fp2.fastq.gz \
        -S $workingdir/bowtie/${i}.sam

        ## Compress the aligned sam files to bam files
        samtools view \
        -bS $workingdir/bowtie/${i}.sam \
        > $workingdir/bowtie/${i}.bam \

        ## Organise mapped reads and index them for faster access during 
        ## downstream processing
        samtools view \
        -@ ${SLURM_CPUS_PER_TASK} \
        -o $workingdir/bowtie/${i}.sorted.bam \
        $workingdir/bowtie/${i}.bam

        samtools index \
        $workingdir/${i}.sorted.bam
        
        ## Run some stats on aligned/mapped reads
        bamtools stats \
        -in $workingdir/bowtie/${i}.sorted.bam \
        > $workingdir/bowtie/${i}.sorted.stats.txt
        ## look at the stats files to quality check the data

done