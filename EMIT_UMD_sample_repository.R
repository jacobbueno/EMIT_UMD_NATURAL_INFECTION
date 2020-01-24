# EMIT_UMD_sample_repository

# Create an authoritative dataset for EMIT UMD natural infection study that includes information on those screened, enrolled, and what the infection was. 
# For some we will not know which pathogen caused illness, but future work looks like it will begin to yield more information about infection pathogen and potential coinfection. 
# It will be helpful to have information also about which samples we have available for future assays, which samples have been sent out to various other labs for analysis, and whether we have received usable data from these other labs or not. 


#---------------------------------

## Set up environment ####

library(tidyverse)

sessionInfo()


## Examine the available analytical datasets and take inventory ####

# Files (located in EMIT/EMIT_Data_Analysis_Jake/EMIT_UMD_Natural_Infection/Curated Data/Analytical Datasets)

# all_cases_gii_samples.csv
# all_cases.csv
# all_screened.csv
# finaldatasetrepeatupdate.csv
# flu_cases_gii_samples_subtype_allignment.csv
# flu_cases_gii_samples.csv
# flu_cases.csv
# PNAS_data_full.csv
# roommates.csv

file_list <- list.files(path = "/Users/jbueno/Box Sync/EMIT/EMIT_Data_Analysis_Jake/EMIT_UMD_Natural_Infection/Curated Data/Analytical Datasets")

files <- lapply(file_list, read.csv)

file_list

names(files) <- file_list

for (i in seq_along(files)) {
  assign(names(files)[i], files[[i]], .GlobalEnv)
}


## For each df, we want to know how many subjects there are, how many g-ii instances are represented 

all_screened <- all_screened.csv %>%
  summarize(n_distinct(subject.id))
all_screened

all_cases <- all_cases.csv %>%
  summarize(n_distinct(subject.id))
all_cases












