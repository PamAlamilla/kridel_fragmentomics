#!/bin/bash
#
#SBATCH -N 1 
#SBATCH -p himem
#SBATCH -c 6
#SBATCH --mem=40000M
#SBATCH -t 0-06:00
#SBATCH -J My-Job #Change to your job name
#SBATCH --array=0-12 #Change to the number of BAM files minus 1

#Description: This script reads through a list of BAM files and applies Samtools view to each one. 
#The output is a new BAM file containing only reads that overlap the genomic regions defined by the BED file.
#Creator: Pamela Alamilla
#Date: 11/06/2024

#Load required module
module load samtools

#Provide path to text file of BAM file names
FILENAMES=/path/to/file.txt #Change this to your text file with the BAM file names

#Provide path to input BAM files
INPUT_PATH=/path/to/directory #Change this to the directory holding your BAM files

#Provide path to directory in which to deposit output BAMs
OUTPUT_PATH=/path/to/directory #Change this to the directory that will hold the output BAMs

#Provide path to BED file of genomic windows
BED_FILE=/path/to/file.bed #Change this to your BED file with the target genomic windows

#Retrieve next BAM file name from provided text file
NAMES=($(cat $FILENAMES))
INPUT_BAM=${NAMES[${SLURM_ARRAY_TASK_ID}]}

#Generate name for output file
OUTPUT_BAM=$(echo "$INPUT_BAM" | sed 's/.\///g' | sed 's/dedup/filtered/g') #Change target pattern to reflect your BAM file names

#Filter BAM file using BED file
samtools view -b -h -M -L $BED_FILE $INPUT_PATH/$INPUT_BAM > $OUTPUT_PATH/$OUTPUT_BAM
