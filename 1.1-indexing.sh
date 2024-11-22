#!/bin/bash

#SBATCH --partition=queue_name      # the requested queue
#SBATCH --nodes=1              # number of nodes to use
#SBATCH --tasks-per-node=1     # for parallel distributed jobs
#SBATCH --cpus-per-task=8      # for multi-threaded jobs
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

module load bowtie2/v2.4.1

## point to the directory containing the reference genome where sequences will be
## mapped

mkdir At_reference_genome/

export refdir=your/working/dir/At_reference_genome

#################################################################################
# Main CMD
#################################################################################

# Retrieving release 59 of the Arabidopsis thaliana reference genome

echo "============================="
echo "Retrieving TAIR10 Release 59"
wget -P $refdir \
        "https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-59/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz"
wget -P $refdir \
        "https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-59/gtf/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.59.gtf.gz"

# Indexing genomes

echo "============================="
echo "RUNNING INDEXING"

## Index the genome for quicker access by bowtie2 during alignment

bowtie2-build $refdir/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz $refdir/Arabidopsis_thaliana.TAIR10.59.gtf.gz


echo "INDEXING COMPLETE"
echo "============================="

#################################################################################
# End
#################################################################################
