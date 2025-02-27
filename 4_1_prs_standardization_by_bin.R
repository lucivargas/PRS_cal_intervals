# Clean environment
rm(list=ls())

# Load libraries
library("rstudioapi") #getActiveDocumentContext
library(tidyverse) # data wrangling
library(modelr) # extend tidyverse
library(BBmisc) # helper functions
library(ggpubr) # data visualization
library("viridis") # add color palettes
library(ggbeeswarm) # data visualization
library(stats) # make auc/roc 
library(epitools) # more epi stats
require(pROC) # evaluate area under the curve
library(DHARMa) # diagnose glm models
library(boot) # bootstrapping for deriving CI
library(BBmisc) # helper functions



########################
# Define file locations 
########################

# Import functions I developed for standardization
source("/home/ldebritovargas@xsede.org/scripts/PRS_standardization_by_bin/4_2_prs_standardization_by_bin_Functions.R")

# Set PRS name
prs_name <- "PGS004603"


# Set working directory
setwd(paste0("/home/ldebritovargas@xsede.org/results/PGS_scores/",prs_name,"/"))

# Phenotype name
#pheno_name <- "Diabetes_SR"

# PRS file
pheno_file <- paste0("/home/ldebritovargas@xsede.org/results/PGS_scores/", prs_name, "/REGARDS_pheno_", prs_name, ".csv")

# Define directory to save results
results_dir <- paste0("/home/ldebritovargas@xsede.org/results/PGS_scores/", prs_name, "/")

########################
# Import & format data
#######################

# Import data
data <- read_csv(pheno_file)
levels(data$AFR_ances_intervals) <- c("EA", "< 50%", "50-60%", 
                                   "60-70%","70-80%", "80-90%", "> 90%")

# Define PRS variables to be standardized
prs_names <- grep(prs_name, names(data), value = T)
prs_names <- prs_names[!grepl("top10", prs_names)]



############################
# Standardize by admx. bin
###########################

# Create variable with admixture intervals 
#data <- create_intervals(data, use.decile=T)

# For each PRS, standardize
for (each_prs in prs_names){

  # Plot before standardization
  plot_intervals(data, prs_value=each_prs, save_plot = T)
  
  # Standardize
  data <- std_intervals(data, prs_value=each_prs)
  
  # Plot after standardization
  plot_intervals(data, prs_value=paste0(each_prs, "_stdxbin"), save_plot = T)
}

# Check scores before and after std by bin
ggplot(data, aes(x=.data[[paste0(prs_name,"_std")]],
                 y=.data[[paste0(prs_name,"_std_stdxbin")]], 
                 color=AFR_ances_prop)) +
  geom_point(alpha=0.3) + viridis::scale_color_viridis(name = "% African Ancestry", discrete = F) +
  geom_abline(lty=2) +
  facet_wrap(vars(Race), nrow=1) + theme_classic()


############################
# Generate at-risk cutoffs
############################

# Flag individuals in top/bottom 25% risk of each PRS, in each population
prs_names <- grep(paste0("(?=.*",prs_name,")(?!.*top)"), names(data), perl=T,  value=T)
for(each_prs in prs_names){
  data <- data %>% 
    group_by(Race) %>% 
    mutate(prs_top25 = .data[[each_prs]] >= 
             quantile(unlist(.data[[each_prs]]), 
                      probs = seq(.25, .75, by = .25), na.rm=T)[3]) %>%
    ungroup() %>%
    mutate(prs_top25 = ifelse(prs_top25, 1, 0)) %>% 
    mutate(across(all_of(c("prs_top25")), as.factor)) %>%
    rename_with(~paste0(each_prs, "_top25"), prs_top25) 
}

# Flag individuals in top/bottom 10% risk of each PRS, in each population
prs_names <- grep(paste0("(?=.*",prs_name,")(?!.*top)"), names(data), perl=T,  value=T)
for(each_prs in prs_names){
  if(paste0(each_prs, "_top10") %in% names(data)){next}
  data <- data %>% 
    group_by(Race) %>% 
    mutate(prs_top10 = .data[[each_prs]] >= 
             quantile(unlist(.data[[each_prs]]), 
                      probs = seq(.1, .9, by = .1), na.rm=T)[9]) %>%
    ungroup() %>%
    mutate(prs_top10 = ifelse(prs_top10, 1, 0)) %>% 
    mutate(across(all_of(c("prs_top10")), as.factor)) %>%
    rename_with(~paste0(each_prs, "_top10"), prs_top10) 
}


# Check that cutoff is correct
ggplot(data, aes(x=.data[[prs_name]], y=id ,color=.data[[paste0(prs_name, "_top10")]])) + 
  geom_point(alpha=0.5) + facet_wrap(vars(Race))

# Get distribution of PRS according to % AFR intervals 
prs_binary_names <- grep("top10", names(data), value = T)
for(each_prs in prs_binary_names){
  prs_distr <- data %>% group_by(AFR_ances_intervals) %>% 
    summarise(count= n(), at_risk = sum(.data[[each_prs]] == 1),
              perc_at_risk = (at_risk/count)*100) %>%
    arrange(match(AFR_ances_intervals, levels(AFR_ances_intervals)))
  print(each_prs)
  print(prs_distr)
}

# Save data

write_csv(data, file=paste0(results_dir, "REGARDS_", prs_name, "_pheno_stdxbin.csv"))

