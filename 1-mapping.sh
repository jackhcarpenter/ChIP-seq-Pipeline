#!/bin/bash

#SBATCH --partition=defq       # the requested queue
#SBATCH --nodes=1              # number of nodes to use
#SBATCH --tasks-per-node=1     # for parallel distributed jobs
#SBATCH --cpus-per-task=8      # for multi-threaded jobs
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

module load bowtie2/v2.4.1
module load samtools/1.10
module load bamtools/v2.5.1

## point to the directory containing the reference genome where sequences will be
## mapped
export refdir=/mnt/scratch/c1831460/ChIP/At_reference_genome

## point to the working directory
export workingdir=/mnt/scratch/c1831460/ChIP

##REMEMBER: set up any directories that the software needs in this script in case 
##it is unable to do so itself

#################################################################################
# Main CMD
#################################################################################

## Index the genome for quicker access by bowtie2 during alignment
#bowtie2-build \
#    $refdir/SPECIES_GENOME.dna.toplevel.fa\
#    $refdir/SPECIES_GENOME.INDEX.fa

## List of sequences to map to indexed reference genome
list=(
	"Col-O-AB_S2" "Col-O-Input_S3" "Col-O-NoAB_S1" "TCPA-AB_S5"\
	"TCPA-Input_S6" "TCPA-NoAB_S4")


## Map forward and reverse reads to the indexed referenced genome
for i in ${list[@]}
do
        
        echo ${i}

        ## Align the sequences to the genome
        bowtie2 \
        --maxins 1000 \
        --fr \
        -p 8 \
	--trim-to 36 \
        -x $refdir/ARABIDOPSIS_THALIANA.GENOME.INDEX.fa \
        -1 $workingdir/RAW_DATA/230628_fastqs/${i}_merge_R1.fastq.gz \
        -2 $workingdir/RAW_DATA/230628_fastqs/${i}_merge_R2.fastq.gz \
        -S $workingdir/bowtie/maxins_1000_trimmed_36/${i}.sam

        ## Compress the aligned sam files to bam files
        samtools view \
        -bS $workingdir/bowtie/maxins_1000_trimmed_36/${i}.sam \
        > $workingdir/bowtie/maxins_1000_trimmed_36/${i}.bam \

        ## Organise mapped reads and index them for faster access during 
        ## downstream processing
        
        samtools sort \
        -@ ${SLURM_CPUS_PER_TASK} \
        -o $workingdir/bowtie/maxins_1000_trimmed_36/${i}.sorted.bam \
        $workingdir/bowtie/maxins_1000_trimmed_36/${i}.bam

        samtools index \
        $workingdir/bowtie/maxins_1000_trimmed_36/${i}.sorted.bam

        ## Run some stats on aligned/mapped reads
        bamtools stats \
        -in $workingdir/bowtie/maxins_1000_trimmed_36/${i}.sorted.bam \
        > $workingdir/bowtie/maxins_1000_trimmed_36/${i}.sorted.stats.txt

        ## look at the stats files to quality check the data

done
