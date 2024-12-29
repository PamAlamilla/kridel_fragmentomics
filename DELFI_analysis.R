setwd("/Users/pamalamilla/Desktop/Fragmentomics_GOOD/DELFI/")

#Combine Bins----
#Lifted and adapted from GitHub: cancer-genomics/delfi_scripts/03.5-combine_bins.r
library(tidyverse)
library(GenomicRanges)
library(multidplyr)
library(readxl)
library(cowplot)

bindir <- "~/Desktop/Fragmentomics_GOOD/DELFI/fragment_profiles_RDS" #Provide directory of RDS files

metadata <- read.table("~/Desktop/Fragmentomics_GOOD/metadata_discovery.txt") #Read in metadata table
metadata <- metadata[2:144,] #Exclude sample LIB_15_0036_T1

metadata_clean <- metadata %>% select(SAMPLE_ID, TYPE, DIAGNOSIS_CLASS)

ids <- metadata_clean %>% select(SAMPLE_ID) %>% unlist() #Create vector of sample IDs

files <- file.path(bindir, paste0(sub("LIB_15_", "LIB-15-", ids),
                                  "_100kb_fragment_profile.RDS")) #Create vector of file paths

bins.list <- lapply(files, readRDS) #Read in RDS files into list

tib.list <- lapply(bins.list, as_tibble) #Convert RDS files to tibbles in list
names(tib.list) <- ids #Name elements of tibbles list

tib.list <- map2(tib.list, names(tib.list), ~ mutate(.x, id = .y)) %>%
  bind_rows() %>% select(id, everything()) #Merge tibbles in list into one tibble
 
tib.list <- tib.list %>% select(-matches("X")) #Exclude columns starting with 'X'

saveRDS(tib.list, "bins_100kbcompartments.rds")

#5mb Bins----
#Lifted and adapted from cancer-genomics/delfi_scripts/04-5mb_bins.r

df.fr <- readRDS("bins_100kbcompartments.rds")

df.fr2 <- inner_join(df.fr, metadata, by=c("id"="SAMPLE_ID"))

armlevels <- c("1p","1q","2p","2q","3p","3q","4p","4q","5p","5q","6p","6q",
               "7p","7q","8p","8q", "9p", "9q","10p","10q","11p","11q","12p",
               "12q","13q","14q","15q","16p","16q","17p","17q","18p","18q",
               "19p", "19q","20p","20q","21q","22q")
df.fr2$arm <- factor(df.fr2$arm, levels=armlevels) #Add column of chromosome arms

## combine adjacent 100kb bins to form 5mb bins. We count starting from
## the telomeric end and remove the bin closest to the centromere if it is
## smaller than 5mb.
df.fr2 <- df.fr2 %>% group_by(id, arm) %>%
  mutate(combine = ifelse(grepl("p", arm), ceiling((1:length(arm))/50),
                          ceiling(rev((1:length(arm))/50) )))

df.fr3 <- df.fr2 %>% group_by(id, seqnames, arm, combine) %>%
  summarize(short2=sum(short),
            long2=sum(long),
            short.corrected2=sum(short.corrected),
            long.corrected2=sum(long.corrected),
            hic.eigen=mean(eigen),
            gc=mean(C.G),
            ratio2=mean(ratio),
            ratio.corrected2=mean(ratio.corrected),
            nfrags2=sum(nfrags),
            nfrags.corrected2=sum(nfrags.corrected),
            domain = median(as.integer(domain)),
            short.var=var(short.corrected),
            long.var=var(long.corrected),
            nfrags.var=var(nfrags.corrected),
            mode_size=unique(mode),
            mean_size=unique(mean),
            median_size=unique(median),
            q25_size=unique(quantile.25),
            q75_size=unique(quantile.75),
            start=start[1],
            end=rev(end)[1],
            binsize = n())

### assign bins
df.fr3 <- inner_join(df.fr3, metadata_clean, by=c("id"="SAMPLE_ID")) #Add TYPE and DIAGNOSIS_CLASS columns

df.fr3 <- df.fr3 %>% filter(binsize==50)

df.fr3 <- df.fr3 %>% group_by(id) %>% mutate(bin = 1:length(id))

saveRDS(df.fr3, "bins_5mbcompartments.rds")


#Summarize Data----
#Lifted and adapted from GitHub: cancer-genomics/delfi_scripts/05-summarize_data.r

df.fr3 <- readRDS("./bins_5mbcompartments.rds")

healthy.median <- df.fr3 %>% group_by(bin) %>% filter(TYPE=="CONTROL") %>% 
  summarize(median.cov=median(nfrags2, na.rm=TRUE),
            median.short=median(short2, na.rm=TRUE),
            median.long=median(long2, na.rm=TRUE),
            median.ratio=median(ratio2, na.rm=TRUE),
            median.corrected.cov=median(nfrags.corrected2, na.rm=TRUE),
            median.corrected.short=median(short.corrected2, na.rm=TRUE),
            median.corrected.long=median(long.corrected2, na.rm=TRUE),
            median.corrected.ratio=median(ratio.corrected2, na.rm=TRUE),
            median.corrected.ratio2=median(short.corrected2/long.corrected2, na.rm=TRUE))

summary.df <- df.fr3 %>% ungroup() %>% group_by(id, TYPE) %>%
  summarize(cov.cor=cor(nfrags2, healthy.median$median.cov, method="pearson", use="complete.obs"),
            short.cor=cor(short2, healthy.median$median.short, method="pearson", use="complete.obs"),
            long.cor=cor(long2, healthy.median$median.long, method="pearson", use="complete.obs"),
            ratio.cor=cor(ratio2, healthy.median$median.ratio, method="pearson", use="complete.obs"),
            cov.corrected.cor=cor(nfrags.corrected2, healthy.median$median.corrected.cov, method="pearson", use="complete.obs"),
            short.corrected.cor=cor(short.corrected2, healthy.median$median.corrected.short, method="pearson", use="complete.obs"),
            long.corrected.cor=cor(long.corrected2, healthy.median$median.corrected.long, method="pearson", use="complete.obs"),
            ratio.corrected.cor=cor(ratio.corrected2, healthy.median$median.corrected.ratio, method="pearson", use="complete.obs"),
            ratio2.corrected.cor=cor(short.corrected2/long.corrected2, healthy.median$median.corrected.ratio2, method="pearson", use="complete.obs"),
            nfrags = sum(nfrags2),
            mode_size=unique(mode_size),
            mean_size=unique(mean_size),
            median_size=unique(median_size),
            q25_size=unique(q25_size),
            q75_size=unique(q75_size),
            hqbases_analyzed = 100*sum(nfrags)*2,
            coverage = hqbases_analyzed/(504*5e6)
  )

summary.df$`TYPE` = relevel(as.factor(summary.df$`TYPE`), "CONTROL")

saveRDS(summary.df, "./summary_tibble.rds")

#Plot Genome-Wide Fragmentation Profiles----
df.fr3 <- readRDS("./bins_5mbcompartments.rds")

df.fr3 <- df.fr3 %>% 
  group_by(id) %>%
  mutate(ratio.centered = scale(ratio.corrected2, scale=FALSE)[,1])

mytheme <- theme_classic(base_size=12) + theme(axis.text.x = element_blank(),
                                               axis.ticks.x=element_blank(),
                                               #strip.text.y = element_blank(), #Uncomment to remove labels along right y-axis
                                               strip.text.x = element_text(face="bold", size=20),
                                               strip.text.y = element_text(face="bold", size=20),
                                               axis.title.x = element_text(face="bold", size=50),
                                               axis.title.y = element_text(face="bold", size=50),
                                               axis.text.y = element_text(size = 20),
                                               plot.title = element_text(size = 15),
                                               legend.position = "none",
                                               panel.grid.major = element_blank(),
                                               panel.grid.minor = element_blank(),
                                               strip.background = element_blank())

tissue <- c(CONTROL = "Healthy (n=?)", LYMPHOMA = "Lymphoma (n=?)")

arm <- df.fr3 %>% group_by(arm) %>%
  summarize(n=n()) %>%
  mutate(arm = as.character(arm))
small.arms <- setNames(c("", "12q", "", "16q", "", "17q", "", "18q","", "19",
                         "", "20", "21", "22"),
                       c("12p", "12q", "16p", "16q", "17p", "17q", "18p", "18q",
                         "19p", "19q", "20p", "20q","21q", "22q"))

arm.labels <- setNames(arm$arm, arm$arm)
arm.labels[names(small.arms)] <- small.arms

df.fr3$DIAGNOSIS_CLASS <- factor(df.fr3$DIAGNOSIS_CLASS)

subtype = c(CONTROL = "CONTROL", FL = "FL", DLBCL = "DLBCL", HL = "HL")

##Controls vs Lymphomas----
gg = ggplot(df.fr3, aes(x=bin, y=ratio.centered, group=id)) +
  geom_line(color="gray50", size=1) +
  labs(x="", y="Fragmentation profile", color="") +
  facet_grid(TYPE~arm, switch="x",space="free_x", scales="free_x", labeller=labeller(type=tissue, arm=arm.labels)) + 
  coord_cartesian(xlim = NULL, ylim=c(-.10,.12), expand = TRUE) +
  mytheme

save_plot("../Figures/DELFI_fragment_profiles_2groups.pdf", gg, ncol=1, nrow=1, base_height=12, base_width=35)

##By Subtype----
gg2 = ggplot(df.fr3, aes(x=bin, y=ratio.centered, group=id)) +
  geom_line(color="gray50", size=1) +
  labs(x="", y="Fragmentation profile", color="") +
  facet_grid(DIAGNOSIS_CLASS~arm, switch="x", space="free_x", scales="free_x", 
             labeller=labeller(DIAGNOSIS_CLASS = subtype, arm=arm.labels)) + 
  coord_cartesian(xlim = NULL, ylim=c(-.10,.12), expand = TRUE) +
  mytheme

save_plot("../Figures/DELFI_fragment_profiles_4groups.pdf", gg2, ncol=1, nrow=1, base_height=24, base_width=35)

#Calculate Profile Correlations----
summary.df.2 = inner_join(summary.df, metadata_clean[,c(1,3)], by=c("id"="SAMPLE_ID")) #Add DIAGNOSIS_CLASS column

control.summary.tib = summary.df %>% filter(TYPE=="CONTROL") %>% select(ratio2.corrected.cor)
lymphoma.summary.tib = summary.df %>% filter(TYPE=="LYMPHOMA") %>% select(ratio2.corrected.cor)
DLBCL.summary.tib = summary.df.2 %>% filter(DIAGNOSIS_CLASS=="DLBCL") %>% select(ratio2.corrected.cor)
FL.summary.tib = summary.df.2 %>% filter(DIAGNOSIS_CLASS=="FL") %>% select(ratio2.corrected.cor)
HL.summary.tib = summary.df.2 %>% filter(DIAGNOSIS_CLASS=="HL") %>% select(ratio2.corrected.cor)

control.summary.ratio2corrected = control.summary.tib[,2] %>% unlist()
lymphoma.summary.ratio2corrected = lymphoma.summary.tib[,2] %>% unlist() 
DLBCL.summary.ratio2corrected = DLBCL.summary.tib[,2] %>% unlist() 
FL.summary.ratio2corrected = FL.summary.tib[,2] %>% unlist() 
HL.summary.ratio2corrected = HL.summary.tib[,2] %>% unlist() 

control.summary.ratio2corrected %>% unlist() %>% median() #0.860458 median correlation
lymphoma.summary.ratio2corrected %>% unlist() %>% median() #0.7616431 median correlation
DLBCL.summary.ratio2corrected %>% unlist() %>% median() #0.7289043 median correlation
FL.summary.ratio2corrected %>% unlist() %>% median() #0.7794597 median correlation
HL.summary.ratio2corrected %>% unlist() %>% median() #0.7665276 median correlation

wilcox.test(control.summary.ratio2corrected, lymphoma.summary.ratio2corrected) #p-value = 1.434e-05
wilcox.test(control.summary.ratio2corrected, DLBCL.summary.ratio2corrected) #p-value = 4.867e-05
wilcox.test(control.summary.ratio2corrected, FL.summary.ratio2corrected) #p-value = 0.004173
wilcox.test(control.summary.ratio2corrected, HL.summary.ratio2corrected) #p-value = 0.0001001


#Plot Corrected Ratios----
##Controls vs Lymphomas----
boxplot(control.summary.ratio2corrected, lymphoma.summary.ratio2corrected, 
        ylim = c(0,1), names = c("Control", "Lymphoma"), 
        main="GC-Corrected Ratios")

##By Subtype----
boxplot(control.summary.ratio2corrected, DLBCL.summary.ratio2corrected, 
        FL.summary.ratio2corrected, HL.summary.ratio2corrected,
        ylim = c(0,1), names = c("Control", "DLBCL", "FL", "HL"),
        main="GC-Corrected Ratios")

DELFI_correctedRatios_4groups
#Dataframe for ML----
df_Victoria = df.fr3 %>% select(id, TYPE, DIAGNOSIS_CLASS, bin, 
                  short.corrected2, long.corrected2) %>%
  pivot_wider( id_cols = c(id, TYPE, DIAGNOSIS_CLASS),
               names_from = bin,
               values_from = c(short.corrected2, long.corrected2))

write_tsv(df_Victoria, "DELFI_fragmentomics.tsv")
