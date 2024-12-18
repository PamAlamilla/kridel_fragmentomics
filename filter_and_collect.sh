#!/bin/bash
#SBATCH -N 1 
#SBATCH -p himem
#SBATCH -c 6
#SBATCH --mem=40000M
#SBATCH -t 0-06:00
#SBATCH -J FilterCollect
#SBATCH --array=0-104 #Change to number of samples minus 1

# Creator: Pamela Alamilla
# Date: 10/10/2024

#Load necessary modules
module load samtools
module load R/4.2.1
module load picard/2.6.0


#Provide path to text file of raw BAM file names
RAW_BAMNAMES=/cluster/projects/kridelgroup/pamela/training_cancers/training_cancers_filenames.txt

#Provide path to input BAM files
RAW_BAMPATH=/cluster/projects/kridelgroup/LIBERATE/KLCS_cfMeDip_2023-2024_completed_samples/All_completed_remove_dup_merge_all_Aug_2024_updated/MEDIPIPE/Aug_2024_total_307_Res/dedup_bam_umi_pe

#Provide path to directory in which to deposit output BAMs
FILTERED_BAMPATH=/cluster/projects/kridelgroup/pamela/training_cancers/DESeq_p0.05_FC1.5/filtered_bams

#Provide path to BED file of genomic windows
BED_FILE=/cluster/projects/kridelgroup/pamela/bed_files/DESeq_p0.05_FC1.5.bed

#Provide paths to directories in which to deposit output txt and pdf files
METRICS_PATH=/cluster/projects/kridelgroup/pamela/training_cancers/DESeq_p0.05_FC1.5/insert_size_metrics
HISTOGRAMS_PATH=/cluster/projects/kridelgroup/pamela/training_cancers/DESeq_p0.05_FC1.5/insert_size_histograms


#Retrieve next BAM file name from provided text file
RAW_NAMES=($(cat $RAW_BAMNAMES))
RAW_BAM=${RAW_NAMES[${SLURM_ARRAY_TASK_ID}]}

#Generate name for output filtered BAM file
FILTERED_BAM=$(echo "$RAW_BAM" | sed 's/.\///g' | sed 's/dedup/filtered/g')

#Filter BAM file using BED file
samtools view -b -h -M -L $BED_FILE $RAW_BAMPATH/$RAW_BAM > $FILTERED_BAMPATH/$FILTERED_BAM

#Generate names for output metrics and histogram files
METRICS=$(echo "$FILTERED_BAM" | sed 's/.\///g' | sed 's/filtered.bam/insert_size_metrics.txt/g')
HISTOGRAM=$(echo "$FILTERED_BAM" | sed 's/.\///g' | sed 's/filtered.bam/insert_size_histogram.pdf/g')

#Collect insert sizes
java -jar $picard_dir/picard.jar CollectInsertSizeMetrics \
	I=$FILTERED_BAMPATH/$FILTERED_BAM \
	O=$METRICS_PATH/$METRICS \
	H=$HISTOGRAMS_PATH/$HISTOGRAM
