# This script merges PRS calculated in a dataset to the phenotype data in that same dataset

# Clean environment
rm(list=ls())

# Load libraries
library(tidyverse)

# Set prs name
prs_name <- "PGS002308"

# Set path to merged weighted PRS scores from step 7
prs_file <- paste0("/scratch/alpine/ldebritovargas@xsede.org/PGS_scoring/",prs_name,"/REGARDS_",prs_name,".csv")


# Set path to phenotype file
pheno_file <- paste0("/home/ldebritovargas@xsede.org/results/PGS_scores/", prs_name, "/REGARDS_pheno.csv")


# Where to save merged file
merged_file <-  paste0("/home/ldebritovargas@xsede.org/results/PGS_scores/",prs_name,"/REGARDS_pheno_",prs_name,".csv")


# Import pheno and PRS
pheno_df <- read_csv(pheno_file) 
prs_df <- read_csv(prs_file) %>% select(-Race) 


# Merge PRS to pheno
if(any(grepl(prs_name, names(pheno_df)))){break}
prs_df <- prs_df %>%
  inner_join(., pheno_df, by="IID") %>%
  relocate(contains(prs_name), .after=contains("afr"))

# Remove individuals where SIRE and %Afr don't match,
# Criteria: Black with %Afr<15%, White that have %afr > 5%
prs_df <- prs_df %>%
  filter(! (Race=="B" & AFR_ances_prop < 0.15),
         ! (Race=="W" & AFR_ances_prop > 0.05))


# Save merged file
write_csv(prs_df, file=merged_file)


# End
