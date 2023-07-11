#!/bin/bash

#SBATCH --partition=QUEUE_NAME       # the requested queue
#SBATCH --nodes=1              # number of nodes to use
#SBATCH --tasks-per-node=1     # for parallel distributed jobs
#SBATCH --cpus-per-task=4      # for multi-threaded jobs
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

## list of samples
list=("sample_input1" "sample_inputn" "sample_ip1" "sample_ipn" "sample_neg1" 
"sample_negn" "control_input1" "control_inputn" "control_ip1" "control_ipn" 
"control_neg1" "control_negn")

## perform fastqc on the raw PE data
for i in list ${list[@]}
do
        echo ${i}

        fastqc $workingdir/${i}_merge_R1.fastq.gz
        fastqc $workingdir/${i}_merge_R2.fastq.gz

done

## summarise the QC data of all reads
multiqc -i "PROJECT_NAME_RAW_SEQUENCES" $workingdir/

## perform fastp to remove low quality reads and adaptors 
for i in list {$list[@]}
do
        echo ${i}

        fastp \
        -i $workingdir/${i}_merge_R1.fastq.gz \
        -I $workingdir/${i}_merge_R2.fastq.gz \
        -o $workingdir/trimmed_reads/${i}_fp1.fastq.gz \
        -O $workingdir/trimmed_reads/${i}_fp2.fastq.gz \
        --detect_adapter_for_pe --trim_poly_g --correction
        
done

## perform fastqc on the trimmed sequences
for i in ${list[@]}
do
        echo ${i}

	    fastqc $workingdir/trimmed_reads/${i}_fp1.fastq.gz
	    fastqc $workingdir/trimmed_reads/${i}_fp2.fastq.gz

done

## summarise the QC data of the filtered and trimmed reads
multiqc -i "PROJECT_NAME" $workingdir/trimmed_reads/