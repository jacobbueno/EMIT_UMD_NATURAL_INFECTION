# Title: Enrollment Condition for Negative Samples
# EMIT_Enroll_Condition_for_Negative_Samples.R
# Author: Jacob Bueno de Mesquita, with material from Jing Yan and Donald Milton

# Summary: I'm moving this script to the git lab repository and also doing some minor cleaning to enabling the clean reproduction of this script.
# Signed: Jacob Bueno de Mesquita
# Date: January 8; February, 2019

# Description: This script reads in the Clinical Database (EMIT_samples.cc.RDS) and the subtyping file (EMIT_subtypes.RDS) to create a summary file that contains a list of the sampling instances with culture data (passage or quantitative focus assay) that shows that there is positive virus, despite a negative first visit NP swab (first visit NP swab is the basis of the EMIT_subtypes.RDS file). 
# The output is: "/Users/jbueno/Box Sync/EMIT/EMIT_Data_Analysis_Jake/EMIT_UMD_Natural_Infection/Curated Data/Cleaned Data/negative subtype sample with positive culture or focus.RDS"
# This output can be used in subsequent analyses or exploration into the dataset as it pertains to counting individuals with evidence of viral infection by both culture and qRT-PCR methods. 

#########################
# Name: enroll condition for negative samples
# Authors: Jing Yan & Donald Milton
# Date: February 3, 2016
# # Purpose:  
# merge "EMIT_subtypes" (df with list of type and subtype from NP swabs) with 
# "samples.cc"  (df with a list of all samples, culture & some clinical data), 
# Then output the list of subjects with negative sample type and alos not enrolled, 
# and negative sample type with positive culture results
#########################

library(dplyr)
library(tidyr)

sessionInfo()

#_____________
# Setup all I/O
# setwd('C:/Users/Jing/Box Sync/EMIT/EMIT_Data_Analysis')
# setwd('/Users/dmilton/Box Sync/0_DKM/Lab/Biodefense/EMIT/EMIT_Data_Analysis')
# setwd('/Volumes/Internal RAID Set 1/Box Sync/0_DKM/Lab/Biodefense/EMIT/EMIT_Data_Analysis')
# setwd("/Users/jbueno/Box Sync/EMIT/EMIT_Data_Analysis_Jake")

# Comment out the working directory to facilitate markdown report compilation

# Out.dir <- "R_output/"
# sink(file = "R_output/enroll condition for negative samples.txt", split = TRUE)
#_____________

samples <- readRDS("/Users/jbueno/Box Sync/EMIT/EMIT_Data_Analysis_Jake/EMIT_UMD_Natural_Infection/Curated Data/Cleaned Data/EMIT_samples.cc.RDS")
flu.types <- readRDS("/Users/jbueno/Box Sync/EMIT/EMIT_Data_Analysis_Jake/EMIT_UMD_Natural_Infection/Curated Data/Cleaned Data/EMIT_subtypes.RDS")
flu.types <- select(flu.types, subject.id, type.inf)

samples.types <- left_join(samples, flu.types, by = "subject.id")

# cat("\n Rows samples:", nrow(samples), " Rows flu.types:", nrow(flu.types), " Rows samples.types:", nrow(samples.types))

subjects.types <- samples.types %>%
  distinct(subject.id)

# print(with(subjects.types,(ftable(addmargins(table(enrolled,type.inf,exclude = c()))))))

neg.samples <- samples.types %>%
  filter(type.inf == 'Negative' | type.inf == 'bad assay')

# cat("\n With a negative sample type (or bad assay), the number of subjects have a screen visit 2 is ",
print(nrow(filter(neg.samples, visit.num == 2)))

# With a negative sample type (or bad assay), the number of subjects have a screen visit 2 is...
nrow(filter(neg.samples, visit.num == 2))

neg.unenroll <- neg.samples %>% 
  filter(g2.run == 0) %>%
  distinct(date.visit, subject.id, .keep_all = TRUE) %>%
  select(subject.id, date.visit, sample.id, sample.type, g2.run, 
                       visit.num, passpos, validp, focus.ct, enrolled)

# cat("\n subjects with negative sample type and are also not enrolled. \n")
print(neg.unenroll) 

enroll.np <- neg.samples %>% 
  filter(g2.run != 0)
neg.pcr.pos.passage <- neg.samples %>% 
  filter(passpos == TRUE & validp == TRUE)
neg.pcr.pos.focus <- neg.samples %>% 
  filter(focus.ct > 0)

pos.culture <- neg.pcr.pos.passage %>%
  full_join(neg.pcr.pos.focus, by = c("subject.id", "type.inf", "date.visit", "sample.id", "sample.type", "visit.num",                                          "passpos", "g2.run", "validp", "focus.ct", "enrolled"))

pos.culture <- pos.culture %>%
  arrange(subject.id) %>% 
  select(subject.id, date.visit, sample.id, sample.type, g2.run, visit.num, passpos, validp, focus.ct, enrolled)

# cat("\n subjects with negative sample type but positive culture results. \n")
print(pos.culture)

saveRDS(pos.culture, "/Users/jbueno/Box Sync/EMIT/EMIT_Data_Analysis_Jake/EMIT_UMD_Natural_Infection/Curated Data/Cleaned Data/negative subtype sample with positive culture or focus.RDS")

# sink()
# closeAllConnections()
