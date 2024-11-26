#Load Packages----
library(ggplot2)
library(reshape2)
library(NMF)
library(stringr)
library(ggsignif)

#Define Functions----
generateSampleNames = function(fileNames) {
  # Takes a vector of file names containing either a LIBERATE or TGL sample ID
  # and returns a sorted vector of sample names.
  
  # Examples: 
  # "LIB-15-0002_T1_insert_size_metrics.txt" becomes "LIB_15_0002_T1"
  # "TGL49_0255_Ct_n_PE_220318_CM_insert_size_metrics.txt" becomes "TGL49_0255_Ct"
  
  sampleNames = str_extract(gsub ("-", "_", fileNames), 
                            "LIB_[:digit:]{2}_[:digit:]{4}_[:alnum:]{2}")
  sampleNames = append(sampleNames, 
                       str_extract(gsub ("-", "_", fileNames), 
                                   "TGL49_[:digit:]{4}_[:alpha:]{2}"))
  sampleNames = sampleNames[!is.na(sampleNames)]
  sampleNames = sort(sampleNames)
  return(sampleNames)
}

collectPercents = function(ListOfSamples, fragment_length) {
  # Takes a list of samples with fragment data and a extracts all the percentages
  # for given fragment_length (missing values are assigned NA).
  
  percents = c()
  for (i in 1:length(ListOfSamples)) {
    percent_vec = as.numeric(subset(ListOfSamples[[i]], insert_size == fragment_length, percentages))
    percents = append(percents,percent_vec)
  }
  return(percents)
}

#____RAW PROFILES____----
# Note: Before starting, download fragment count files from h4h cluster to local computer. Sort into 
# two directories, one for control files and another for lymphoma files. Fragment count files 
# are formatted like so: <path_to_MEDIPIPE>/MEDIPIPE/Res/fragment_size/Sample_ID_insert_size_metrics.txt

#Read in fragment count tables----
##Controls----
setwd("./controls_raw") #Move to directory with fragment count files for controls

# Make vector of file names to be read in
controls_fileNames = sort(list.files(pattern = "*_insert_size_metrics.txt"))

# Make vector of sample names
controls_sampleNames = generateSampleNames(controls_fileNames)
length(controls_sampleNames) == length(controls_fileNames) #Expect TRUE

# Make list of sample data & assign names to samples
controls_raw_list = lapply(controls_fileNames, FUN=read.table, skip=10, header=T)
names(controls_raw_list) = controls_sampleNames

#Add column of fraction percentages to each sample in list
for (i in 1:length(controls_raw_list)) {
  totalCount = sum(controls_raw_list[[i]][2])
  percentages_vector = (controls_raw_list[[i]][[2]]/totalCount) *100
  controls_raw_list[[i]]["percentages"] <- percentages_vector
}
