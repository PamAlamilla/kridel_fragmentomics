#!/bin/bash
#
#SBATCH -N 1 
#SBATCH -p himem
#SBATCH -c 6
#SBATCH --mem=40000M
#SBATCH -t 0-06:00
#SBATCH -J My_Job #Change to your job name
#SBATCH --array=0-32 #Change to the number of BAM files minus 1

#Description: This script reads through a list of BAM files and applies the CollectInsertSizeMetrics tool (from Picard) to each one. 
#The output is two files, a text file containing insert size metrics and a PDF file with the associated histogram. 
#Author: Pamela Alamilla
#Date: 05/06/2024

#Load necessary modules
module load R/4.2.1
module load picard/2.6.0

#Provide path to text file of BAM file names
FILENAMES=/path/to/file.txt #Change this to your file with the BAM file names

#Provide path to BAM files
BAM_PATH=/path/to/directory #Change this to the directory holding your BAM files

#Provide paths to directories in which to deposit output files
METRICS_PATH=/path/to/directory #Change this to the directory that will hold the output text files
HISTOGRAMS_PATH=/path/to/directory #Change this to the directory that will hold the output PDFs

#Retrieve next BAM file name from provided text file
NAMES=($(cat $FILENAMES))
INPUT_BAM=${NAMES[${SLURM_ARRAY_TASK_ID}]}

#Generate names for output files
METRICS=$(echo "$INPUT_BAM" | sed 's/.\///g' | sed 's/filtered.bam/insert_size_metrics.txt/g') #Change target pattern to reflect your BAM file names
HISTOGRAM=$(echo "$INPUT_BAM" | sed 's/.\///g' | sed 's/filtered.bam/insert_size_histogram.pdf/g') #Change target pattern to reflect your BAM file names

#Collect insert sizes
java -jar $picard_dir/picard.jar CollectInsertSizeMetrics \
	I=$BAM_PATH/$INPUT_BAM \
	O=$METRICS_PATH/$METRICS \
	H=$HISTOGRAMS_PATH/$HISTOGRAM
