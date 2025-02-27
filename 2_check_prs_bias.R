# Clean environment
rm(list=ls())

# Load libraries
library(tidyverse)
library(ggridges)
library(ggpubr)

# Set up
setwd("/home/ldebritovargas@xsede.org/results/PGS_scores/")

########################
# Define file locations 
########################

# Set PRS name and outcome of interest
prs_name <- "PGS002308"

# Set paths to clean data set with PRS of interest
pheno_file <- paste0("/home/ldebritovargas@xsede.org/results/PGS_scores/",prs_name,"/REGARDS_pheno_",prs_name,".csv")

# Where to save plots
results_folder <- paste0("/home/ldebritovargas@xsede.org/results/PGS_scores/", prs_name, "/img/")
if(!dir.exists(results_folder)){dir.create(results_folder)}

###############
## Import data
###############

# Import pheno+PRS
x <- read_csv(pheno_file) 
levels(x$AFR_ances_intervals) <- c("0%", "< 50%", "50-60%", 
                                        "60-70%","70-80%", "80-90%", "> 90%")

# Rename PRS
x <- x %>%
  rename_with(~ gsub(prs_name, "grs.wt", .x, fixed = TRUE)) 


################
## Check data
################

# Subset by SIRE 
x_afr <- x[x$Race == "B",]
x_eur <- x[x$Race == "W",]
# & then by sex
x_afr_F <- x_afr[x_afr$Gender == "F",]
x_afr_M <- x_afr[x_afr$Gender == "M",]
x_eur_F <- x_eur[x_eur$Gender == "F",]
x_eur_M <- x_eur[x_eur$Gender == "M",]


# Check correlation of PRS and % AA ancestry in AAs
ggplot(x_afr, aes(x=grs.wt_std, y=AFR_ances_prop)) + 
  geom_point(alpha=0.3, color="lightblue3") + 
  geom_smooth(method=lm, level=0.95, color="blue4", alpha=0.7) + stat_cor() + 
  labs(x=prs_name, y="% of African ancestry") + theme_classic()
ggsave(paste0(results_folder, "img/REGARDS_GRS-",prs_name,"_cor_ancestry.png"), height = 5, width = 7)
summary(lm(grs.wt_std ~ AFR_ances_prop, data=x_afr))


# Check 'at-risk' numbers per group - BEFORE std_by_bin
# check size of each bin
table(x_afr$AFR_ances_intervals)
# check distribution of high-risk in each bin
table(x_afr[,c("grs.wt_std_top10","AFR_ances_intervals")])
# check percentage of high-risk in each bin
(table(x_afr[,c("grs.wt_std_top10","AFR_ances_intervals")])[2,]/table(x_afr$AFR_ances_intervals))*100
# check % of high-risk in whole cohort
table(x_afr$grs.wt_top10)[2]/sum(table(x_afr$grs.wt_top10))


# End
