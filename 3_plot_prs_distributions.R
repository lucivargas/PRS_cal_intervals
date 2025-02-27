# Clean environment
rm(list=ls())

# Load libraries
library(tidyverse)
library(ggridges)
library(ggpubr)

# Set PRS name 
prs_name <- "PGS002308"

# Set up directory
setwd(paste0("/home/ldebritovargas@xsede.org/results/PGS_scores/",prs_name,"/"))

########################
# Define file locations 
########################

# Set paths to data set with weighted PRS scores from PRS evaluation step 1 (merging)
pheno_file <- paste0("/home/ldebritovargas@xsede.org/results/PGS_scores/",prs_name,"/REGARDS_pheno_",prs_name,".csv")

# Where to save plots
results_folder <- paste0("/home/ldebritovargas@xsede.org/results/PGS_scores/", prs_name, "/img/")


###############
## Import data
###############

# Create results dir if one doesnt exist
if(!dir.exists(results_folder)){dir.create(results_folder)}

# Import data
x <- read_csv(pheno_file)

# Rename PRS
x <- x %>%
  rename_with(~ gsub(prs_name, "grs.wt", .x, fixed = TRUE)) 

# Reorder ancestry intervals
x <- x %>%
  mutate(AFR_ances_intervals = factor(AFR_ances_intervals, levels = c("EA", "< 50%", "50-60%", "60-70%", "70-80%", "80-90%", "> 90%")))

################
## Create plots
################

# Plot PRS in REGARDS by SIRE
# grs.wt
ggplot(x, aes(x=grs.wt, color=Race)) +
  scale_color_discrete(labels = c("African-American", "European-American")) +
  labs(y = "Density",  x = expression("PRS"[T2D])) +
  geom_density() +
  theme_classic()
ggsave(paste0(results_folder, "REGARDS_",prs_name,"_by_SIRE.png"), height = 5, width = 7)
# grs.wt_std
ggplot(x, aes(x=grs.wt_std, color=Race)) +
  scale_color_discrete(labels = c("African-American", "European-American")) +
  labs(y = "Density",  x = expression("PRS"[T2D])) +
  geom_density() +
  theme_classic()
ggsave(paste0(results_folder, "REGARDS_",prs_name,"_std_by_SIRE.png"), height = 5, width = 7)

# Plot PRS in REGARDS by admixture classes
ggplot(x, aes(x=grs.wt, color=AFR_ances_intervals)) +
  viridis::scale_color_viridis(name = "% African Ancestry", discrete = TRUE) +
  #scale_color_discrete(labels = c("African-American", "European-American")) +
  labs(y = "Density",  x = expression("PRS"[T2D])) +
  geom_density(size=1.2) +
  theme_classic() 
ggsave(paste0(results_folder, "REGARDS_",prs_name,"_by_ancestry.png"), height = 5, width = 8.5)

# Grouped Scatter plot with marginal density plots
ggscatterhist(
  x, x = "grs.wt", y = "AFR_ances_prop", color = "Race", 
  margin.plot = "histogram", size = 3, alpha = 0.2,
  palette = c("red", "#00AFBB"),
  main.plot.size = 1, margin.plot.size = 2,
  margin.params = list(fill = "Race", color = "black", size = 0.2, n = 2),
  ylab = "% African Ancestry",  xlab = "T2D PRS",
  legend = "none") 
ggsave(paste0(results_folder, "REGARDS_",prs_name,"_vs_ancestry.png"), height = 7, width = 10)

# End
