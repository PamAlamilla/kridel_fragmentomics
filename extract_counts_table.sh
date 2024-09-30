#!/bin/bash

#Description: This simple script takes an insert_size_metrics text file output by MEDIPIPE and extracts the table of fragment counts (headings included).
#The output is a new text file with the tab-separated table.
#Creator: Pamela Alamilla
#Date: 06/06/2024

for FILE in path/to/directory/*_insert_size_metrics.txt #Change this to the directory holding your metrics files
do
	OUTPUT_FILE=$(echo "$FILE" | sed 's/.\///g' | sed 's/insert_size_metrics/fragment_counts/g') #Note the new file names have different endings.
	cat $FILE | tail -n +11 > $OUTPUT_FILE
done
