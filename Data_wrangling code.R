# Load necessary library
library(dplyr)
library(tidyr)

# Read the CSV file
antibiotics <- read.csv("2024_05_28 atlas_antibiotics.csv")

# List of genotype columns to keep
genotype_columns <- c("AMPC", "SHV", "TEM", "CTXM1", "CTXM2", "CTXM825", "CTXM9", "VEB", "PER", 
                      "GES", "ACC", "CMY1MOX", "CMY11", "DHA", "FOX", "ACTMIR", "KPC", "OXA", 
                      "NDM", "IMP", "VIM", "SPM", "GIM")

# Remove columns that do not end with '_I', but keep the genotype columns
antibiotics <- antibiotics %>%
  select(Isolate.Id:Phenotype, all_of(genotype_columns), ends_with("_I"))

names(antibiotics)

# Convert from wide to long format for the antibiotic columns
antibiotics_long <- antibiotics %>%
  pivot_longer(cols = ends_with("_I"),
               names_to = "Antibiotic",
               values_to = "Resistance") %>%
  drop_na(Resistance)  # Exclude rows with missing Resistance values


