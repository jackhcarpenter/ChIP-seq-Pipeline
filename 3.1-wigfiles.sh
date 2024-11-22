#!/bin/bash

#SBATCH --partition=queue_name       # the requested queue
#SBATCH --nodes=1              # number of nodes to use
#SBATCH --tasks-per-node=1     # for parallel distributed jobs
#SBATCH --cpus-per-task=8      # for multi-threaded jobs
#SBATCH --mem-per-cpu=2G      # in megabytes, unless unit explicitly stated
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

## point to the directory containing the reference genome where sequences will be
## mapped
export refdir=your/working/dir/At_reference_genome

## point to the working directory
export workingdir=/your/working/dir

mkdir beds
mkdir wigs

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

for file in ${files[@]}
do
        echo ${file} "= running..."

## Use samtools to create an index in fasta format so as betools can access quickly
    samtools faidx $refdir/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa

## Bin (make discrete) map sequences to create histograms
        bedtools genomecov \
        -ibam $workingdir/bowtie/${file}.sorted.bam \
        -bg \
        -g $redir/Arabidopsis_thaliana.release59.TAIR10.dna.toplevel.fa.fai > \
        $workingdir/beds/${file}.bedgraph

        echo ${file} "genome covarage = complete"

## Generate a coverage track to view in Integrated Genome Viewer
        bamCoverage \
        -b $workingdir/bowtie/${file}.sorted.bam \
        -o $workingdir/wigs/${file}.bw

        echo ${file} "bam covarage = complete"

done
#################################################################################
# End
#################################################################################