#!/bin/bash

#SBATCH --partition=queue_name       # the requested queue
#SBATCH --nodes=1              # number of nodes to use
#SBATCH --tasks-per-node=1     #
#SBATCH --cpus-per-task=4      #
#SBATCH --mem-per-cpu=4G       # in megabytes, unless unit explicitly stated
#SBATCH --error=logs/%J.err         # redirect stderr to this file
#SBATCH --output=logs/%J.out        # redirect stdout to this file
#SBATCH --mail-user=your.email@host  # email
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

module load picard/2.22.2
module load bamtools/v2.5.1
module load samtools/1.10

# point to the directory containing the reference genome

export refdir=your/working/dir/At_referene_genome
echo "Reference directory =" $refdir

# define the working directory

export workingdir=your/working/dir
echo "Working directory =" $workingdir

# define the export directory

exportdir=/mnt/scratch/c1831460/ChIP/markdup
echo "Export directory =" $exportdir

##REMEMBER: set up any directories that the software needs in this script in case
##it is unable to do so itself

#################################################################################
# Main CMD
#################################################################################

# Loop variables

TF=("true" \
        "false")

declare -a files

for file in $workingdir/bowtie/*sorted.bam
do
        echo ${file}
        files+=("$(basename ${file::-11})")

done

# Loops

for var in ${TF[@]}
do

        for file in ${files[@]}
        do
                if ${var} == "true"
                then
                        echo ${file} "remove duplicate = running"
                        dup="rmdup"
                else
                        echo ${file} "mark duplicate = running"
                        dup="markdup"
                fi

                ## Remove the duplicated reads
                java -jar $PICARD MarkDuplicates \
                        I=$workingdir/bowtie/${file}.sorted.bam \
                        O=$exportdir/${file}_${dup}.bam \
                        M=$exportdir/${file}_metrics.${dup}.txt \
                        REMOVE_DUPLICATES=${var} \
                        VALIDATION_STRINGENCY=SILENT

                bamtools stats \
                        -in $exportdir/${file}_${dup}.bam \
                        > $exportdir/${file}_${dup}.dupstats.txt

                if ${var} == "true"
                then
                        echo ${file} "remove duplicate = complete"
                else
                        echo ${file} "mark duplicate = complete"
                fi

                echo "============================="

                ## Now look at the files to see if it is better to keep or remove duplicates

        done
done
#################################################################################
# End
#################################################################################