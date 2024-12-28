setwd("/Users/pamalamilla/Desktop/Fragmentomics_GOOD/")

#Load Packages----
library(ggplot2)
library(reshape2)
library(NMF)
library(stringr)
library(ggsignif)
library(ggpubr)

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
  # Takes a list of samples with fragment data and a extracts all the fractions
  # for given fragment_length (missing values are assigned NA).
  percents = c()
  for (i in 1:length(ListOfSamples)) {
    percent_vec = as.numeric(subset(ListOfSamples[[i]], insert_size == fragment_length, percentages))
    percents = append(percents,percent_vec)
  }
  return(percents)
}

calculatePPS = function(matrix){
  # Takes a matrix of fragmentation profiles and calculates the per-patient score 
  # of each sample, per method published by Vessies et al (20##). Matrix must 
  # contain samples as rows and 900 rows (i.e. fragment sizes 1bp to 900bp). 
  # Table of fragment scores must also have been read in.
  
  patient_scores = c()
  for (i in 1:nrow(matrix)) {
    fragments_sample = sample(1:900, 1000000, replace=T, prob=matrix[i,])
    FS = mean(Vessies_FS[fragments_sample,1])
    patient_scores = append(patient_scores, FS)
  }
  return(patient_scores)
}

#Read in Vessies FS Table----
Vessies_FS = read.table("./Vessies_Fragmentation_Score.csv")

#Read in Metadata Table----
metadata = read.table("./metadata_discovery.txt")

#Define Colours for Plots----
lymphomaRED_transparent = rgb(red=1,green=0,blue=0, alpha=0.3)
controlBLUE_transparent = rgb(red=0,green=0.749,blue=1, alpha=0.3)

#********************************----
#RAW PROFILES----
#********************************----
# Note: Before starting, download fragment count files from h4h cluster to local computer. Sort into 
# two directories, one called "training_controls_raw" and another called "training_lymphomas_raw". 
# Fragment count files are formatted like so:
# <path_to_MEDIPIPE>/MEDIPIPE/Res/fragment_size/Sample_ID_insert_size_metrics.txt

#Read in Fragmentation Profiles----
##Controls----
setwd("./training_controls_raw") #Move to directory with raw fragmentation profiles for controls

controls_fileNames = sort(list.files(pattern = "*_insert_size_metrics.txt")) #Make vector of file names to be read in
controls_sampleNames = generateSampleNames(controls_fileNames) #Make vector of sample names
length(controls_sampleNames) == length(controls_fileNames) #Expect TRUE
controls_raw_list = lapply(controls_fileNames, FUN=read.table, skip=10, header=T) #Make list of sample data
names(controls_raw_list) = controls_sampleNames #Assign names to samples

#Add column of fraction percentages to each sample in list
for (i in 1:length(controls_raw_list)) {
  totalCount = sum(controls_raw_list[[i]][2])
  percentages_vector = (controls_raw_list[[i]][[2]]/totalCount) *100
  controls_raw_list[[i]]["percentages"] <- percentages_vector
}

##Lymphomas----
setwd("../training_lymphomas_raw") #Move to directory with raw fragmentation profiles for lymphomas

lymphomas_fileNames = sort(list.files(pattern = "*_insert_size_metrics.txt")) #Make vector of file names to be read in
lymphomas_sampleNames = generateSampleNames(lymphomas_fileNames) #Make vector of sample names
length(lymphomas_sampleNames) == length(lymphomas_fileNames) #Expect TRUE
lymphomas_raw_list = lapply(lymphomas_fileNames, FUN=read.table, skip=10, header=T) #Make list of sample data
names(lymphomas_raw_list) = lymphomas_sampleNames #Assign names to samples

#Add column of fraction percentages to each sample in list
for (i in 1:length(lymphomas_raw_list)) {
  totalCount = sum(lymphomas_raw_list[[i]][2])
  percentages_vector = (lymphomas_raw_list[[i]][[2]]/totalCount) *100
  lymphomas_raw_list[[i]]["percentages"] <- percentages_vector
}


#Make Matrices----
##Controls----
controls_raw_mtrx = matrix(nrow = length(controls_raw_list), ncol = 900)
rownames(controls_raw_mtrx) = controls_sampleNames #Name rows by samples
colnames(controls_raw_mtrx) = 1:ncol(controls_raw_mtrx)

for (fragment_length in 1:ncol(controls_raw_mtrx)) {
  controls_raw_mtrx[,fragment_length] <- collectPercents(controls_raw_list, fragment_length)
}

controls_raw_mtrx[is.na(controls_raw_mtrx)] = 0 #Set NA values to 0

##Lymphomas----
lymphomas_raw_mtrx = matrix(nrow = length(lymphomas_raw_list), ncol = 900)
rownames(lymphomas_raw_mtrx) = lymphomas_sampleNames #Name rows by samples
colnames(lymphomas_raw_mtrx) = 1:ncol(lymphomas_raw_mtrx)

for (fragment_length in 1:ncol(lymphomas_raw_mtrx)) {
  lymphomas_raw_mtrx[,fragment_length] <- collectPercents(lymphomas_raw_list, fragment_length)
}

lymphomas_raw_mtrx[is.na(lymphomas_raw_mtrx)] = 0 #Set NA values to 0

###DLBCL----
setwd("..")
DLBCL_sample_names = unlist(read.table("training_DLBCL_49.txt"))
DLBCL_sample_indices = which(rownames(lymphomas_raw_mtrx) %in% DLBCL_sample_names)
DLBCL_raw_mtrx = lymphomas_raw_mtrx[which(rownames(lymphomas_raw_mtrx) 
                                          %in% DLBCL_sample_names),]

###FL----
FL_sample_names = unlist(read.table("training_FL_30.txt"))
FL_sample_indices = which(rownames(lymphomas_raw_mtrx) %in% FL_sample_names)
FL_raw_mtrx = lymphomas_raw_mtrx[which(rownames(lymphomas_raw_mtrx) 
                                       %in% FL_sample_names),]

###HL----
HL_sample_names = unlist(read.table("training_HL_31.txt"))
HL_sample_indices = which(rownames(lymphomas_raw_mtrx) %in% HL_sample_names)
HL_raw_mtrx = lymphomas_raw_mtrx[which(rownames(lymphomas_raw_mtrx) 
                                       %in% HL_sample_names),]


#Calculate Mean + SD----
controls_raw_mean = apply(controls_raw_mtrx, 2, mean)
controls_raw_sd = apply(controls_raw_mtrx, 2, sd)

lymphomas_raw_mean = apply(lymphomas_raw_mtrx, 2, mean)
lymphomas_raw_sd = apply(lymphomas_raw_mtrx, 2, sd)

#Generate Plots----
##Controls Only----
ggplot(data=melt(controls_raw_mtrx[,30:250]), 
       aes(x=Var2, y=value, colour=Var1)) +
  geom_line() +
  theme_bw() +
  ggtitle("Controls Raw Fragmentation Profiles") +
  xlab("Fragment Size (bp)") +
  ylab("Frequency (%)") +
  expand_limits(x = c(30, 250), y = c(0,4.0)) +
  scale_x_continuous(breaks = seq(from=30, to=250, by = 20)) +
  scale_y_continuous(breaks = seq(from=0, to=4.0, by = 0.5)) +
  theme(legend.position="none") +
  scale_color_manual(values = rep(controlBLUE_transparent, nrow(controls_raw_mtrx)))

##Lymphomas Only----
ggplot(data=melt(lymphomas_raw_mtrx[,30:250]), 
       aes(x=Var2, y=value, colour=Var1)) +
  geom_line() +
  theme_bw() +
  ggtitle("Lymphomas Raw Fragmentation Profiles") +
  xlab("Fragment Size (bp)") +
  ylab("Frequency (%)") +
  expand_limits(x = c(30, 250), y = c(0,4.0)) +
  scale_x_continuous(breaks = seq(from=30, to=250, by = 20)) +
  scale_y_continuous(breaks = seq(from=0, to=4.0, by = 0.5)) +
  theme(legend.position="none") +
  scale_color_manual(values = rep(lymphomaRED_transparent, nrow(lymphomas_raw_mtrx)))

##Mean Fragmentation Profiles----
ggplot() +
  geom_ribbon(aes(x=30:250, y = controls_raw_mean[30:250], 
                  ymin = (controls_raw_mean - controls_raw_sd)[30:250], 
                  ymax = (controls_raw_mean + controls_raw_sd)[30:250],), 
              fill = controlBLUE_transparent) +
  geom_ribbon(aes(x=30:250, y = lymphomas_raw_mean[30:250], 
                  ymin = (lymphomas_raw_mean - lymphomas_raw_sd)[30:250], 
                  ymax = (lymphomas_raw_mean + lymphomas_raw_sd)[30:250]), 
              fill = lymphomaRED_transparent) +
  geom_line(aes(x=30:250, y = controls_raw_mean[30:250], linetype="Control"), colour = "blue") +
  geom_line(aes(x=30:250, y = lymphomas_raw_mean[30:250], linetype="Lymphoma"), colour = "red") +
  theme_bw() +
  ggtitle("Mean Raw Fragmentation Profiles") +
  xlab("Fragment Size (bp)") +
  ylab("Frequency (%)") +
  expand_limits(x = c(30, 250), y = c(0,3.5)) +
  scale_x_continuous(breaks = seq(from=30, to=250, by = 20)) +
  scale_y_continuous(breaks = seq(from=0, to=3.5, by = 0.5)) +
  scale_linetype_manual(name='Legend',
                        breaks=c("Control", "Lymphoma"),
                        values=c("Control"=1, "Lymphoma"=1)) + 
  theme(legend.title = element_blank())

#Vessies Per-Patient Score----
##Calculate Scores----
set.seed(1)
controls_raw_PPS = calculatePPS(controls_raw_mtrx)
lymphomas_raw_PPS = calculatePPS(lymphomas_raw_mtrx)

##Calculate Cutoff----
raw_PPS_cutoff = mean(controls_raw_PPS) + sd(controls_raw_PPS)*2

##Sensitivity + Specificity----
sum(lymphomas_raw_PPS > raw_PPS_cutoff)/length(lymphomas_raw_PPS)*100 #32/110 or 29.09%
sum(controls_raw_PPS <= raw_PPS_cutoff)/length(controls_raw_PPS) *100 #32/33 or 94.12%

##Plot----
###Controls vs Lymphomas----
raw_PPS_matrix = melt(cbind.data.frame("SAMPLE"=append(rep("Control",nrow(controls_raw_mtrx)), 
                                                       rep("Lymphoma",nrow(lymphomas_raw_mtrx))),
                             "FS"=append(controls_raw_PPS, lymphomas_raw_PPS)), id.vars="SAMPLE")

ggplot(raw_PPS_matrix, aes(x=SAMPLE,y=value, fill=SAMPLE)) +
  geom_boxplot() +
  geom_point(size=1) +
  theme_bw() +
  ylim(-0.75,1)+
  ggtitle("Vessies Per-Patient Score (Raw)") +
  xlab("Sample Type") +
  ylab("Score") +
  geom_signif(comparisons = list(c("Control", "Lymphoma"))) +
  geom_hline(yintercept=raw_PPS_cutoff, linetype="dashed") + 
  scale_fill_manual(values=c(controlBLUE_transparent, lymphomaRED_transparent))+
  theme(legend.position="none")

###By Subtype----
ggplot() +
  geom_boxplot(aes(x="Control", y=controls_raw_PPS), fill=controlBLUE_transparent) +
  geom_point(aes(x="Control", y=controls_raw_PPS), col="black", size=1) +
  
  geom_boxplot(aes(x="DLBCL", y=lymphomas_raw_PPS[DLBCL_sample_indices]),fill="#48F5D1") +
  geom_point(aes(x="DLBCL", y=lymphomas_raw_PPS[DLBCL_sample_indices]),col="black", size=1) +
  
  geom_boxplot(aes(x="FL", y=lymphomas_raw_PPS[FL_sample_indices]),fill="#51F360") +
  geom_point(aes(x="FL", y=lymphomas_raw_PPS[FL_sample_indices]),col="black", size=1) +
  
  geom_boxplot(aes(x="HL", y=lymphomas_raw_PPS[HL_sample_indices]),fill="#FFC60C", ) +
  geom_point(aes(x="HL", y=lymphomas_raw_PPS[HL_sample_indices]),col="black", size=1) +
  
  theme_bw() +
  ylim(-0.75,1)+
  ggtitle("Per-Patient Scores by Subtype (Raw)") +
  xlab("Sample Type") +
  ylab("Score") +
  geom_hline(yintercept=raw_PPS_cutoff, linetype="dashed")


#NMF Rank 2----
set.seed(2)
NMF2_raw_results = nmf(rbind(controls_raw_mtrx[,30:250],
                             lymphomas_raw_mtrx[,30:250]), 2)

##Scale Signatures----
NMF2_raw_Sig1_scaled = c()
for (i in 1:ncol(NMF2_raw_results@fit@H)) {
  h_matrix = NMF2_raw_results@fit@H
  sample_value = h_matrix[1,i]/sum(h_matrix[1,])
  NMF2_raw_Sig1_scaled = append(NMF2_raw_Sig1_scaled, sample_value)
}

NMF2_raw_Sig2_scaled = c()
for (i in 1:ncol(NMF2_raw_results@fit@H)) {
  h_matrix = NMF2_raw_results@fit@H
  sample_value = h_matrix[2,i]/sum(h_matrix[2,])
  NMF2_raw_Sig2_scaled = append(NMF2_raw_Sig2_scaled, sample_value)
}

##Plot Signatures----
ggplot() +
  geom_line(aes(x=30:250,
                y=NMF2_raw_Sig1_scaled*100,
                linetype="Control Signature")) +
  geom_line(aes(x=30:250,
                y=NMF2_raw_Sig2_scaled*100,
                linetype="Lymphoma Signature")) +
  theme_bw() +
  ggtitle("NMF Rank 2 Signatures (Raw)")+
  expand_limits(x = c(30, 250), y = c(0,3.5)) +
  scale_x_continuous(breaks = seq(from=30, to=250, by = 20)) +
  scale_y_continuous(breaks = seq(from=0, to=3.5, by = 0.5)) +
  xlab("Fragment Size (bp)") +
  ylab("Frequency (%)") +
  scale_linetype_manual(name='Legend',
                        breaks=c("Control Signature", "Lymphoma Signature"),
                        values=c("Control Signature"=1, "Lymphoma Signature"=5)) + 
  theme(legend.title = element_blank())

#********************************----
#HYPERMETHYLATED PROFILES----
#********************************----
#Read in Fragmentation Profiles----
##Controls----
setwd("./training_controls_hypermethylated") #Move to directory with hypermethylated fragmentation profiles for controls

#NOTE: Sample LIB-15-0165_T1 is missing from this analysis, hence the c(1:11,13:34) part
controls_hypermethylated_list = lapply(controls_fileNames[c(1:11,13:34)], FUN=read.table, skip=10, header=T) #Make list of sample data
names(controls_hypermethylated_list) = controls_sampleNames[c(1:11,13:34)] #Assign names to samples

#Add column of percentages
for (i in 1:length(controls_hypermethylated_list)) {
  totalCount = sum(controls_hypermethylated_list[[i]][2])
  percentages_vector = (controls_hypermethylated_list[[i]][[2]]/totalCount) *100
  controls_hypermethylated_list[[i]]["percentages"] <- percentages_vector
}

##Lymphomas----
setwd("../training_lymphomas_hypermethylated/") #Move to directory with hypermethylated fragmentation profiles for lymphomas

available_files = c(seq(from = 1, to = 43, by = 2),45:48,54:60,62:80,83:86,93,101,108) #Remove this part once all the hypermethylated profiles have been generated
lymphomas_hypermethylated_list = lapply(lymphomas_fileNames[available_files], FUN=read.table, skip=10, header=T) #Make list of sample data
names(lymphomas_hypermethylated_list) = lymphomas_sampleNames[available_files] #Assign names to samples

#Add column of percentages
for (i in 1:length(lymphomas_hypermethylated_list)) {
  totalCount = sum(lymphomas_hypermethylated_list[[i]][2])
  percentages_vector = (lymphomas_hypermethylated_list[[i]][[2]]/totalCount) *100
  lymphomas_hypermethylated_list[[i]]["percentages"] <- percentages_vector
}


#Make Matrices----
##Controls----
controls_hypermethylated_mtrx = matrix(nrow = length(controls_hypermethylated_list), ncol = 900)
rownames(controls_hypermethylated_mtrx) = controls_sampleNames[c(1:11,13:34)] #Name rows by samples
colnames(controls_hypermethylated_mtrx) = 1:900

for (fragment_length in 1:ncol(controls_hypermethylated_mtrx)) {
  controls_hypermethylated_mtrx[,fragment_length] <- collectPercents(controls_hypermethylated_list, fragment_length)
}

controls_hypermethylated_mtrx[is.na(controls_hypermethylated_mtrx)] = 0 #Set NA values to 0

##Lymphomas----
lymphomas_hypermethylated_mtrx = matrix(nrow = length(lymphomas_hypermethylated_list), ncol = 900)
rownames(lymphomas_hypermethylated_mtrx) = lymphomas_sampleNames[available_files] #Name rows by samples
colnames(lymphomas_hypermethylated_mtrx) = 1:900

for (fragment_length in 1:ncol(lymphomas_hypermethylated_mtrx)) {
  lymphomas_hypermethylated_mtrx[,fragment_length] <- collectPercents(lymphomas_hypermethylated_list, fragment_length)
}

lymphomas_hypermethylated_mtrx[is.na(lymphomas_hypermethylated_mtrx)] = 0 #Set NA values to 0


#Calculate Mean + SD----
controls_hypermethylated_mean = apply(controls_hypermethylated_mtrx, 2, mean)
controls_hypermethylated_sd = apply(controls_hypermethylated_mtrx, 2, sd)

lymphomas_hypermethylated_mean = apply(lymphomas_hypermethylated_mtrx, 2, mean)
lymphomas_hypermethylated_sd = apply(lymphomas_hypermethylated_mtrx, 2, sd)

#Generate Plots----
##Controls Only----
ggplot(data=melt(controls_hypermethylated_mtrx[,30:250]), 
       aes(x=Var2, y=value, colour=Var1)) +
  geom_line() +
  theme_bw() +
  ggtitle("Controls Hypermethylated Fragmentation Profiles") +
  xlab("Fragment Size (bp)") +
  ylab("Frequency (%)") +
  expand_limits(x = c(30, 250), y = c(0,4.0)) +
  scale_x_continuous(breaks = seq(from=30, to=250, by = 20)) +
  scale_y_continuous(breaks = seq(from=0, to=4.0, by = 0.5)) +
  theme(legend.position="none") +
  scale_color_manual(values = rep(controlBLUE_transparent, nrow(controls_hypermethylated_mtrx)))

##Lymphomas Only----
ggplot(data=melt(lymphomas_hypermethylated_mtrx[,30:250]), 
       aes(x=Var2, y=value, colour=Var1)) +
  geom_line() +
  theme_bw() +
  ggtitle("Lymphomas Hypermethylated Fragmentation Profiles") +
  xlab("Fragment Size (bp)") +
  ylab("Frequency (%)") +
  expand_limits(x = c(30, 250), y = c(0,4.0)) +
  scale_x_continuous(breaks = seq(from=30, to=250, by = 20)) +
  scale_y_continuous(breaks = seq(from=0, to=4.0, by = 0.5)) +
  theme(legend.position="none") +
  scale_color_manual(values = rep(lymphomaRED_transparent, nrow(lymphomas_hypermethylated_mtrx)))

##Mean Fragmentation Profiles----
ggplot() +
  geom_ribbon(aes(x=30:250, y = controls_hypermethylated_mean[30:250], 
                  ymin = (controls_hypermethylated_mean - controls_hypermethylated_sd)[30:250], 
                  ymax = (controls_hypermethylated_mean + controls_hypermethylated_sd)[30:250],), 
              fill = controlBLUE_transparent) +
  geom_ribbon(aes(x=30:250, y = lymphomas_hypermethylated_mean[30:250], 
                  ymin = (lymphomas_hypermethylated_mean - lymphomas_hypermethylated_sd)[30:250], 
                  ymax = (lymphomas_hypermethylated_mean + lymphomas_hypermethylated_sd)[30:250]), 
              fill = lymphomaRED_transparent) +
  geom_line(aes(x=30:250, y = controls_hypermethylated_mean[30:250], linetype="Control"), colour = "blue") +
  geom_line(aes(x=30:250, y = lymphomas_hypermethylated_mean[30:250], linetype="Lymphoma"), colour = "red") +
  theme_bw() +
  ggtitle("Mean Hypermethylated Fragmentation Profiles") +
  xlab("Fragment Size (bp)") +
  ylab("Frequency (%)") +
  expand_limits(x = c(30, 250), y = c(0,3.5)) +
  scale_x_continuous(breaks = seq(from=30, to=250, by = 20)) +
  scale_y_continuous(breaks = seq(from=0, to=3.5, by = 0.5)) +
  scale_linetype_manual(name='Legend',
                        breaks=c("Control", "Lymphoma"),
                        values=c("Control"=1, "Lymphoma"=1)) + 
  theme(legend.title = element_blank())


#Vessies Per-Patient Score----
##Calculate Score----
set.seed(1)
controls_hypermethylated_PPS = calculatePPS(controls_hypermethylated_mtrx)
lymphomas_hypermethylated_PPS = calculatePPS(lymphomas_hypermethylated_mtrx)

##Calculate Cutoff----
hypermethylated_PPS_cutoff = mean(controls_hypermethylated_PPS) + sd(controls_hypermethylated_PPS)*2 #-0.3094095

##Sensitivity + Specificity----
sum(lymphomas_hypermethylated_PPS > hypermethylated_PPS_cutoff)/length(lymphomas_hypermethylated_PPS)*100
sum(controls_hypermethylated_PPS <= hypermethylated_PPS_cutoff)/length(controls_hypermethylated_PPS)*100

##Plot----
###Controls vs Lymphomas----
hypermethylated_PPS_matrix = melt(cbind.data.frame("SAMPLE"=append(rep("Control",33), rep("Lymphoma",59)),
                                                "FS"=append(controls_hypermethylated_PPS, lymphomas_hypermethylated_PPS)), id.vars="SAMPLE")

ggplot(hypermethylated_PPS_matrix, aes(x=SAMPLE,y=value, fill=SAMPLE)) +
  geom_boxplot() +
  geom_point(size=1) +
  theme_bw() +
  ylim(-0.75,1)+
  ggtitle("Vessies Per-Patient Score (Hypermethylated)") +
  xlab("Sample Type") +
  ylab("Score") +
  geom_signif(comparisons = list(c("Control", "Lymphoma"))) +
  geom_hline(yintercept=hypermethylated_PPS_cutoff, linetype="dashed") + 
  scale_fill_manual(values=c(controlBLUE_transparent, lymphomaRED_transparent))+
  theme(legend.position="none")

###By Subtype----
ggplot() +
  geom_boxplot(aes(x="Control", y=controls_hypermethylated_PPS), fill=controlBLUE_transparent) +
  geom_point(aes(x="Control", y=controls_hypermethylated_PPS), col="black", size=1) +
  
  geom_boxplot(aes(x="DLBCL", y=lymphomas_hypermethylated_PPS[DLBCL_sample_indices]),fill="#48F5D1") +
  geom_point(aes(x="DLBCL", y=lymphomas_hypermethylated_PPS[DLBCL_sample_indices]),col="black", size=1) +
  
  geom_boxplot(aes(x="FL", y=lymphomas_hypermethylated_PPS[FL_sample_indices]),fill="#51F360") +
  geom_point(aes(x="FL", y=lymphomas_hypermethylated_PPS[FL_sample_indices]),col="black", size=1) +
  
  geom_boxplot(aes(x="HL", y=lymphomas_hypermethylated_PPS[HL_sample_indices]),fill="#FFC60C", ) +
  geom_point(aes(x="HL", y=lymphomas_hypermethylated_PPS[HL_sample_indices]),col="black", size=1) +

  theme_bw() +
  ylim(-0.75,1)+
  ggtitle("Per-Patient Scores by Subtype (Hypermethylated)") +
  xlab("Sample Type") +
  ylab("Score") +
  geom_hline(yintercept=hypermethylated_PPS_cutoff, linetype="dashed")

#Mann Whitney U ???????????????????----
for (i in 1:ncol(controls_hypermethylated_mtrx)) {
  p_value = wilcox.test(controls_hypermethylated_mtrx[,i], lymphomas_hypermethylated_mtrx[,i])$p.value
  if (p_value != "NA" & p_value != "NaN" & p_value < 0.05) {
      print(i)
  }
}


#NMF Rank 2----
#Recall: NMF will not work on rows of zeros, so these must be excluded.
set.seed(2)
NMF2_hypermethylated_results = nmf(rbind(controls_hypermethylated_mtrx[c(1:22,24:33),30:250],
                               lymphomas_hypermethylated_mtrx[,30:250]), 2)

##Scale Signatures----
NMF2_hypermethylated_Sig1_scaled = c()
for (i in 1:ncol(NMF2_hypermethylated_results@fit@H)) {
  h_matrix = NMF2_hypermethylated_results@fit@H
  sample_value = h_matrix[1,i]/sum(h_matrix[1,])
  NMF2_hypermethylated_Sig1_scaled = append(NMF2_hypermethylated_Sig1_scaled, sample_value)
}

NMF2_hypermethylated_Sig2_scaled = c()
for (i in 1:ncol(NMF2_hypermethylated_results@fit@H)) {
  h_matrix = NMF2_hypermethylated_results@fit@H
  sample_value = h_matrix[2,i]/sum(h_matrix[2,])
  NMF2_hypermethylated_Sig2_scaled = append(NMF2_hypermethylated_Sig2_scaled, sample_value)
}

##Plot Signatures----
ggplot() +
  geom_line(aes(x=30:250,
                y=NMF2_hypermethylated_Sig2_scaled*100,
                linetype="Control Signature")) +
  geom_line(aes(x=30:250,
                y=NMF2_hypermethylated_Sig1_scaled*100,
                linetype="Lymphoma Signature")) +
  theme_bw() +
  ggtitle("NMF Rank 2 Signatures (Hypermethylated)")+
  expand_limits(x = c(30, 250), y = c(0,3.5)) +
  scale_x_continuous(breaks = seq(from=30, to=250, by = 20)) +
  scale_y_continuous(breaks = seq(from=0, to=3.5, by = 0.5)) +
  xlab("Fragment Size (bp)") +
  ylab("Frequency (%)") +
  scale_linetype_manual(name='Legend',
                        breaks=c("Control Signature", "Lymphoma Signature"),
                        values=c("Control Signature"=1, "Lymphoma Signature"=5)) + 
  theme(legend.title = element_blank())

#********************************----
#MY SCORE----
#********************************----
#Raw NMF Signatures----
##Calculate Fragment Scores----
NMF2_raw_differences = NMF2_raw_Sig2_scaled*100 - NMF2_raw_Sig1_scaled*100

###Plot Fragment Scores----
ggplot() +
  geom_line(aes(x=30:250,y=NMF2_raw_differences))+
  theme_bw() +
  ggtitle("My Raw Fragment Scores")+
  scale_x_continuous(breaks = seq(from=30, to=250, by = 20)) +
  scale_y_continuous(breaks = seq(from=-2, to=2, by = 0.5)) +
  xlab("Fragment Length (bp)") +
  ylab("Frequency Difference (%)") +
  geom_hline(yintercept=0, color = "red") 

##Generate Per-Patient Scores----
controls_raw_scores = c()
for (sample in 1:nrow(controls_raw_mtrx)) {
  patient_score = sum(controls_raw_mtrx[sample,30:150] * (NMF2_raw_differences[1:121]))
  controls_raw_scores = append(controls_raw_scores, patient_score)
}

lymphomas_raw_scores = c()
for (sample in 1:nrow(lymphomas_raw_mtrx)) {
  patient_score = sum(lymphomas_raw_mtrx[sample,30:150] * (NMF2_raw_differences[1:121]))
  lymphomas_raw_scores = append(lymphomas_raw_scores, patient_score)
}

##Calculate Cutoff----
raw_score_cutoff = mean(controls_raw_scores) + sd(controls_raw_scores)*2

##Sensitivity + Specificity----
sum(controls_raw_scores <= raw_score_cutoff)/length(controls_raw_scores) #32/34 or 94.12% specificity
sum(lymphomas_raw_scores > raw_score_cutoff)/length(lymphomas_raw_scores) #34/59 or 30.1% specificity

##Plot----
ggplot() +
  ylim(0,70) +
  geom_boxplot(aes(x="Control", y=controls_raw_scores), fill=controlBLUE_transparent) +
  geom_point(aes(x="Control", y=controls_raw_scores), col="black") +
  geom_boxplot(aes(x="Lymphoma", y=lymphomas_raw_scores),fill=lymphomaRED_transparent) +
  geom_point(aes(x="Lymphoma", y=lymphomas_raw_scores),col="black") +
  theme_bw() +
  ggtitle("My Per-Patient Scores (Raw)") +
  xlab("Sample Type") +
  ylab("Score") +
  geom_hline(yintercept=raw_score_cutoff, linetype="dashed")

#Hypermethylated NMF Signatures----
##Calculate Fragment Scores----
NMF2_hypermethylated_differences = NMF2_hypermethylated_Sig1_scaled*100 - NMF2_hypermethylated_Sig2_scaled*100

###Plot Fragment Scores----
ggplot() +
  geom_line(aes(x=30:250,
                y=NMF2_hypermethylated_differences))+
  theme_bw() +
  ggtitle("My Hypermethylated Fragment Scores")+
  scale_x_continuous(breaks = seq(from=30, to=250, by = 20)) +
  scale_y_continuous(breaks = seq(from=-2, to=2, by = 0.5)) +
  xlab("Fragment Length (bp)") +
  ylab("Frequency Difference (%)") +
  geom_hline(yintercept=0, color = "red") 

##Generate Per-Patient Scores----
controls_hypermethylated_scores = c()
for (sample in 1:nrow(controls_hypermethylated_mtrx)) {
  patient_score = sum(controls_hypermethylated_mtrx[sample,30:150] * (NMF2_hypermethylated_differences[1:121]))
  controls_hypermethylated_scores = append(controls_hypermethylated_scores, patient_score)
}

lymphomas_hypermethylated_scores = c()
for (sample in 1:nrow(lymphomas_hypermethylated_mtrx)) {
  patient_score = sum(lymphomas_hypermethylated_mtrx[sample,30:150] * (NMF2_hypermethylated_differences[1:121]))
  lymphomas_hypermethylated_scores = append(lymphomas_hypermethylated_scores, patient_score)
}

##Calculate Cutoff----
hypermethylated_score_cutoff = mean(controls_hypermethylated_scores) + sd(controls_hypermethylated_scores)*2

##Sensitivity + Specificity----
sum(controls_hypermethylated_scores <= hypermethylated_score_cutoff)/length(controls_hypermethylated_scores) #33/33 or 100% specificity
sum(lymphomas_hypermethylated_scores > hypermethylated_score_cutoff)/length(lymphomas_hypermethylated_scores) #32/59 or 54.24% specificity

##Plot----
ggplot() +
  ylim(0,70) +
  geom_boxplot(aes(x="Control", y=controls_hypermethylated_scores), fill=controlBLUE_transparent) +
  geom_point(aes(x="Control", y=controls_hypermethylated_scores), col="black") +
  geom_boxplot(aes(x="Lymphoma", y=lymphomas_hypermethylated_scores),fill=lymphomaRED_transparent) +
  geom_point(aes(x="Lymphoma", y=lymphomas_hypermethylated_scores),col="black") +
  theme_bw() +
  ggtitle("My Per-Patient Scores (Hypermethylated)") +
  xlab("Sample Type") +
  ylab("Score") +
  geom_hline(yintercept=hypermethylated_score_cutoff, linetype="dashed")

#Compare My PPS to Vessies PPS----
ggplot() +
  geom_line(aes(x=30:250,
                y=(NMF2_hypermethylated_Sig1_scaled*100 - NMF2_hypermethylated_Sig2_scaled*100)))+
  geom_line(aes(x=30:250,
                y=(NMF2_raw_Sig2_scaled*100 - NMF2_raw_Sig1_scaled*100)), col="violet")+
  geom_line(aes(x=30:250,
                y=Vessies_FS[30:250,1]), col="orange")+
  theme_bw() +
  ggtitle("Fragment Scores")+
  scale_x_continuous(breaks = seq(from=30, to=250, by = 20)) +
  scale_y_continuous(breaks = seq(from=-2, to=2, by = 0.5)) +
  xlab("Fragment Length (bp)") +
  ylab("Difference in Relative Frequency") +
  geom_hline(yintercept=0, color = "red")


#********************************----
#FIGURE FOR PUBLICATION----
#********************************----
plotsTheme = theme(
  legend.position = "none",
  legend.title = element_blank(),
  axis.title = element_text(size=10))

raw_profiles = ggplot() +
  geom_ribbon(aes(x=30:250, y = controls_raw_mean[30:250], 
                  ymin = (controls_raw_mean - controls_raw_sd)[30:250], 
                  ymax = (controls_raw_mean + controls_raw_sd)[30:250],), 
              fill = lymphomaRED_transparent) +
  geom_ribbon(aes(x=30:250, y = lymphomas_raw_mean[30:250], 
                  ymin = (lymphomas_raw_mean - lymphomas_raw_sd)[30:250], 
                  ymax = (lymphomas_raw_mean + lymphomas_raw_sd)[30:250]), 
              fill = controlBLUE_transparent) +
  geom_line(aes(x=30:250, y = controls_raw_mean[30:250], linetype="Control"), colour = "red") +
  geom_line(aes(x=30:250, y = lymphomas_raw_mean[30:250], linetype="Lymphoma"), colour = "deepskyblue") +
  theme_bw() +
  xlab("Fragment Size (bp)") +
  ylab("Frequency (%)") +
  expand_limits(x = c(30, 250), y = c(0,3.5)) +
  scale_x_continuous(breaks = seq(from=30, to=250, by = 20)) +
  scale_y_continuous(breaks = seq(from=0, to=3.5, by = 0.5)) +
  scale_linetype_manual(name='',
                        breaks=c("Control", "Lymphoma"),
                        values=c("Control"=1, "Lymphoma"=1)) +
  plotsTheme

hypermethylated_profiles = ggplot() +
  geom_ribbon(aes(x=30:250, y = controls_hypermethylated_mean[30:250], 
                  ymin = (controls_hypermethylated_mean - controls_hypermethylated_sd)[30:250], 
                  ymax = (controls_hypermethylated_mean + controls_hypermethylated_sd)[30:250],), 
              fill = lymphomaRED_transparent) +
  geom_ribbon(aes(x=30:250, y = lymphomas_raw_mean[30:250], 
                  ymin = (lymphomas_hypermethylated_mean - lymphomas_hypermethylated_sd)[30:250], 
                  ymax = (lymphomas_hypermethylated_mean + lymphomas_hypermethylated_sd)[30:250]), 
              fill = controlBLUE_transparent) +
  geom_line(aes(x=30:250, y = controls_hypermethylated_mean[30:250], linetype="Control"), colour = "red") +
  geom_line(aes(x=30:250, y = lymphomas_hypermethylated_mean[30:250], linetype="Lymphoma"), colour = "deepskyblue") +
  theme_bw() +
  xlab("Fragment Size (bp)") +
  ylab("Frequency (%)") +
  expand_limits(x = c(30, 250), y = c(0,3.5)) +
  scale_x_continuous(breaks = seq(from=30, to=250, by = 20)) +
  scale_y_continuous(breaks = seq(from=0, to=3.5, by = 0.5)) +
  scale_linetype_manual(name='Legend',
                        breaks=c("Control", "Lymphoma"),
                        values=c("Control"=1, "Lymphoma"=1)) +
  plotsTheme

raw_NMF = ggplot() +
  plotsTheme +
  geom_line(aes(x=30:250,
                y=NMF2_raw_Sig1_scaled*100,
                linetype="Control Signature")) +
  geom_line(aes(x=30:250,
                y=NMF2_raw_Sig2_scaled*100,
                linetype="Lymphoma Signature")) +
  theme_bw() +
  expand_limits(x = c(30, 250), y = c(0,3.5)) +
  scale_x_continuous(breaks = seq(from=30, to=250, by = 20)) +
  scale_y_continuous(breaks = seq(from=0, to=3.5, by = 0.5)) +
  xlab("Fragment Size (bp)") +
  ylab("Frequency (%)") +
  scale_linetype_manual(name="",
                        breaks=c("Control Signature", "Lymphoma Signature"),
                        values=c("Control Signature"=1, "Lymphoma Signature"=5)) +
  plotsTheme

hypermethylated_NMF = ggplot() +
  plotsTheme +
  geom_line(aes(x=30:250,
                y=NMF2_hypermethylated_Sig2_scaled*100,
                linetype="Control Signature")) +
  geom_line(aes(x=30:250,
                y=NMF2_hypermethylated_Sig1_scaled*100,
                linetype="Lymphoma Signature")) +
  theme_bw() +
  expand_limits(x = c(30, 250), y = c(0,3.5)) +
  scale_x_continuous(breaks = seq(from=30, to=250, by = 20)) +
  scale_y_continuous(breaks = seq(from=0, to=3.5, by = 0.5)) +
  xlab("Fragment Size (bp)") +
  ylab("Frequency (%)") +
  scale_linetype_manual(name="",
                        breaks=c("Control Signature", "Lymphoma Signature"),
                        values=c("Control Signature"=1, "Lymphoma Signature"=5)) + 
  plotsTheme
figure <- ggarrange(raw_profiles, raw_NMF, hypermethylated_profiles, hypermethylated_NMF,
                    labels = c("A", "B", "C", "D"),
                    ncol = 2, nrow = 2)
figure
