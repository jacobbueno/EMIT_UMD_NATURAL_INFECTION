# EMIT_UMD_Natural_Infection_Analysis.R
# Program Objective: Run analysis on EMIT Natural Infection datasets
# Author: Jacob Bueno de Mesquita using material from Jing Yan and Don Milton
# Date: December 14, 2018 -- February 2019, August 2019
# Summary: While most of the analysis for the EMIT_UMD study was done in SAS (especially that as part of the PNAS manuscript), here are some others analyses that we will include as well. Here, we estimate the mean and variance for the ratio of ffu/copy number in fine aerosols from EMIT using a mixed model with random effect of person and fixed effect of day since onset of symptoms. 

#### Load required packages and set working directory ####

library(tidyverse)
library(RcppRoll)
library(readxl)
library(knitr)
library(data.table)
library(lubridate)
library(arsenal)
library(lme4)

# setwd("/Users/jbueno/Box Sync/EMIT/EMIT_Data_Analysis_Jake/EMIT_UMD_Natural_Infection")
# Commented out the setwd() command in order to properly run the markdown report compilation.

sessionInfo() # for reproducibility

#### Estimate the mean and variance for the ratio of ffu/copy number in fine aerosols from EMIT ####
# Use a mixed model with random effect of person and fixed effect of day since onset of symptoms. 
# First have to cut the right df

PNAS_data_full <- read.csv("/Users/jbueno/Box Sync/EMIT/EMIT_Data_Analysis_Jake/EMIT_UMD_Natural_Infection/Curated Data/Analytical Datasets/PNAS_data_full.csv")
PNAS_data_full$date.visit <- as.Date(PNAS_data_full$date.visit)

# Find each sample where we have PCR and focus assay

focus_subjects <- PNAS_data_full %>%
  filter(!is.na(focus.ct)) %>%
  distinct(subject.id, date.visit, .keep_all = TRUE) %>%
  distinct(subject.id)

focus <- PNAS_data_full %>%
  filter(!is.na(focus.ct))

focus_pos_subjects <- focus %>%
  filter(focus.ct > 0) %>%
  distinct(subject.id, date.visit, .keep_all = TRUE) %>%
  distinct(subject.id)

focus_pos <- focus %>%
  filter(focus.ct > 0)

focus_pos_fine_subjects <- focus_pos %>%
  filter(sample.type == "GII condensate NO mask") %>%
  distinct(subject.id, date.visit) %>%
  distinct(subject.id)
print(nrow(focus_pos_fine_subjects))

focus_pos_person_visits <- focus_pos %>%
  distinct(subject.id, date.visit)
print(nrow(focus_pos_person_visits))

focus_pos_fine <- focus_pos %>%
  filter(sample.type == "GII condensate NO mask")

# (Could also find the number of fine positive where focus was negative (make sure this uses the set where focus was taken))

# Find the number focus positive where fine was negative (both replicates)
focus_pos_fine_ever_pos_person_visits <- focus_pos_fine %>%
  filter(!is.na(final.copies)) # %>%
# distinct(subject.id, date.visit)
# print(nrow(focus_pos_fine_ever_pos_person_visits))

focus_pos_fine_ever_pos_subjects <- focus_pos_fine %>%
  filter(!is.na(final.copies)) %>%
  distinct(subject.id)
print(nrow(focus_pos_fine_ever_pos_subjects))

focus_pos_fine_neg <- focus_pos_fine %>%
  anti_join(focus_pos_fine_ever_pos_subjects) %>%
  distinct(subject.id, date.visit)
print(nrow(focus_pos_fine_neg))

focus_pos_fine_pos <- focus_pos_fine %>%
  anti_join(focus_pos_fine_neg)

focus_pos_fine_pos_person_visits <- focus_pos_fine %>%
  anti_join(focus_pos_fine_neg) %>%
  distinct(subject.id, date.visit)
print(nrow(focus_pos_fine_pos_person_visits))

# Find the set of sample instances where focus was positive and fine was 1 positive and 1 negative
# Use the LOD/sqrt(2) to impute a value for these nondetects. 
# Looks like there are 6 instances where there are NA for final.copies (5 fluA and 1 fluB)
# focus_pos_fine_pos$final.copies[is.na(focus_pos_fine_pos$final.copies)] <- (500/(2)^(.5))

# Now we are changing the plan and will exclude pcr results if one of the replicates was a nondetect
one_or_both_replicates_negative <- focus_pos_fine %>%
  filter(is.na(final.copies)) %>%
  distinct(subject.id, date.visit)

focus_pos_fine_pos <- focus_pos_fine %>%
  anti_join(one_or_both_replicates_negative)

focus_pos_fine_pos_averaged_replicates <- focus_pos_fine_pos %>%
  group_by(subject.id, date.visit) %>%
  mutate(avg.final.copies = mean(final.copies)) %>%
  mutate(avg.focus.ct = mean(focus.ct)) %>%
  mutate(copy_number_ffu_ratio = avg.final.copies / avg.focus.ct) %>%
  mutate(log10_avg.final.copies = log10(avg.final.copies)) %>%
  mutate(log10_avg.focus.ct = log10(avg.focus.ct)) %>%
  mutate(log10_copy_number_log10_ffu_ratio = log10_avg.final.copies / log10_avg.focus.ct) %>%
  mutate(log10_ratio = log10(copy_number_ffu_ratio)) %>%
  distinct(subject.id, date.visit, .keep_all = TRUE)

## Estimate the mean and variance for the ratio of ffu/copy number in fine aerosols.
## Use a mixed model with random effect of person and fixed effect of day since onset of symptoms. 

# Random intercept (to take by-subject variability into account) with fixed mean model 
# (assumes the effect of dpo to be the same across all subjects)
model_1 <- lmer(copy_number_ffu_ratio ~ dpo + (1|subject.id), data = focus_pos_fine_pos_averaged_replicates)
summary(model_1)
coef(model_1)

summary(focus_pos_fine_pos_averaged_replicates$copy_number_ffu_ratio)
# This is quite skewed

model_1_log_scale <- lmer(log10_copy_number_log10_ffu_ratio ~ dpo + (1|subject.id), data = focus_pos_fine_pos_averaged_replicates)
summary(model_1_log_scale)

summary(focus_pos_fine_pos_averaged_replicates$log10_copy_number_log10_ffu_ratio)
# This is less skewed but the team decided to go with a different plan - see below.

# But if we want the fixed effect of dpo to not be continuous, then we will convert this dpo variable to character
# focus_pos_fine_pos_averaged_replicates$dpo <- ordered(focus_pos_fine_pos_averaged_replicates$dpo)
# Could use dpo as ordinal - as done above, but decided to go with using dpo as a character variable for now. 

focus_pos_fine_pos_averaged_replicates$dpo <- as.character(focus_pos_fine_pos_averaged_replicates$dpo)

# Also, in fitting this second model, we want to use the log10(ratio) because the model was somewhat skewed before
model_2 <- lmer(log10_ratio ~ dpo + (1|subject.id), data = focus_pos_fine_pos_averaged_replicates)
summary(model_2)

summary(focus_pos_fine_pos_averaged_replicates$log10_ratio)

table(focus_pos_fine_pos_averaged_replicates$dpo)

# writing this out so that Charles Ma can work on it related to the power calculation for new RO1 grants due Februrary 5, 2019
# write.csv(focus_pos_fine_pos_averaged_replicates, "/Users/jbueno/Desktop/focus_pos_fine_pos.csv")



## For Gene Tan JCVI Sequencing ####

## What were the ct values of the influenza B virus samples (among positive replicates)

# Flu B NPS

flub_nps <- PNAS_data_full %>%
  filter(type == "B") %>%
  filter(sample.type == "Nasopharyngeal swab") %>%
  filter(!is.na(Ct)) %>%
  summarize(n_fluB_nps = n(), mean_ct = mean(Ct)) 
print(flub_nps)

# Flu B fine aerosol

flub_fine <- PNAS_data_full %>%
  filter(type == "B") %>%
  filter(sample.type == "GII condensate NO mask") %>%
  filter(!is.na(Ct)) %>%
  summarize(n_fluB_fine = n(), mean_ct = mean(Ct)) 
print(flub_fine)


# Flu B coarse aerosol

flub_coarse <- PNAS_data_full %>%
  filter(type == "B") %>%
  filter(sample.type == "Impactor 5 um NO mask") %>%
  filter(!is.na(Ct)) %>%
  summarize(n_fluB_coarse = n(), mean_ct = mean(Ct)) 
print(flub_coarse)


# Flu A NPS

flua_nps <- PNAS_data_full %>%
  filter(type == "A") %>%
  filter(sample.type == "Nasopharyngeal swab") %>%
  filter(!is.na(Ct)) %>%
  summarize(n_fluA_nps = n(), mean_ct = mean(Ct)) 
print(flua_nps)


# Flu A fine aerosol

flua_fine <- PNAS_data_full %>%
  filter(type == "A") %>%
  filter(sample.type == "GII condensate NO mask") %>%
  filter(!is.na(Ct)) %>%
  summarize(n_fluA_fine = n(), mean_ct = mean(Ct)) 
print(flua_fine)


# Flu A coarse aerosol

flua_coarse <- PNAS_data_full %>%
  filter(type == "A") %>%
  filter(sample.type == "Impactor 5 um NO mask") %>%
  filter(!is.na(Ct)) %>%
  summarize(n_fluA_coarse = n(), mean_ct = mean(Ct)) 
print(flua_coarse)


## Checking JCVI sequencing metadata with Todd Treagen - August 2019

## Read in the batch data in the files from JCVI (these link up with the sequence results)

UMD_batch1 <- read_xls("/Users/jbueno/Box Sync/EMIT/EMIT_Data_Analysis_Jake/EMIT_UMD_Natural_Infection/Sequence_Data/JCVI/20190306_UMDA_batch1_final_stats.xls", sheet = "UMDA_batch1")

UMD_batch1_draft_status <- read_xls("/Users/jbueno/Box Sync/EMIT/EMIT_Data_Analysis_Jake/EMIT_UMD_Natural_Infection/Sequence_Data/JCVI/20190306_UMDA_batch1_final_stats.xls", sheet = "Draft Status", col_names = FALSE)

UMD_batch2 <- read_xls("/Users/jbueno/Box Sync/EMIT/EMIT_Data_Analysis_Jake/EMIT_UMD_Natural_Infection/Sequence_Data/JCVI/20190411_UMDA_batch2_final_stats.xls", sheet = "UMDA_batch2")

UMD_batch2_draft_status <- read_xls("/Users/jbueno/Box Sync/EMIT/EMIT_Data_Analysis_Jake/EMIT_UMD_Natural_Infection/Sequence_Data/JCVI/20190411_UMDA_batch2_final_stats.xls", sheet = "Draft Status", col_names = FALSE)


## Deal with the formatting in the "Draft Status" sheets -- fill in the JCVI Bac ID for each place where it is missing (in the UMD_batch1_draft_status and UMD_batch2_draft_status dfs)

## For batch1 and batch2, set the column names, remove the blank rows, then fill in the JCVI Bac IDs where needed

# batch1...

UMD_batch1_draft_status_cols <- UMD_batch1_draft_status %>%
  rename(JCVI_Bac_ID = X__1,
         segment = X__2,
         note = X__3) %>%
  filter(!is.na(segment))

for (row in 2:length(UMD_batch1_draft_status_cols$JCVI_Bac_ID)) { # 2 to not affect column names
  if(is.na(UMD_batch1_draft_status_cols$JCVI_Bac_ID[row])) {  
    UMD_batch1_draft_status_cols$JCVI_Bac_ID[row] = UMD_batch1_draft_status_cols$JCVI_Bac_ID[row - 1] 
  }
}


# batch2...


UMD_batch2_draft_status_cols <- UMD_batch2_draft_status %>%
  rename(JCVI_Bac_ID = X__1,
         segment = X__2,
         note = X__3) %>%
  filter(!is.na(segment))

for (row in 2:length(UMD_batch2_draft_status_cols$JCVI_Bac_ID)) { # 2 to not affect column names
  if(is.na(UMD_batch2_draft_status_cols$JCVI_Bac_ID[row])) {  
    UMD_batch2_draft_status_cols$JCVI_Bac_ID[row] = UMD_batch2_draft_status_cols$JCVI_Bac_ID[row - 1] 
  }
}


## Merge batch data with the draft status for batches 1 and 2 and join these batches together by binding rows

UMD_batch1_with_status <- UMD_batch1 %>%
  left_join(UMD_batch1_draft_status_cols, by = c("JCVI Bac ID" = "JCVI_Bac_ID"))

UMD_batch2_with_status <- UMD_batch2 %>%
  left_join(UMD_batch2_draft_status_cols, by = c("JCVI Bac ID" = "JCVI_Bac_ID"))

UMD_batch_1_and_2_with_status <- UMD_batch1_with_status %>%
  bind_rows(UMD_batch2_with_status)


## Read in the sample metadata produced at UMD

metadata_sample_codes_fluA_raw <- read.csv("/Users/jbueno/Box Sync/EMIT/EMIT_Data_Analysis_Jake/EMIT_UMD_Natural_Infection/Sequence_Data/JCVI/Metadata/20180923_NIGSP_UMDA_00001-00239.csv")

metadata_sample_codes_fluB_raw <- read.csv("/Users/jbueno/Box Sync/EMIT/EMIT_Data_Analysis_Jake/EMIT_UMD_Natural_Infection/Sequence_Data/JCVI/Metadata/20180923_NIGSP_UMDB_00001-00089.csv")


## Clean and merge the fluA and fluB metadata

metadata_sample_codes_fluA <- metadata_sample_codes_fluA_raw %>%
  filter(Attribute == "") %>%
  select(Blinded.Number, Tube.Label, Virus.Name, Serotype, Collection.Date, Comments, Special.Notes, Date.Sent.to.JCVI, BioProject.ID) %>%
  filter(Blinded.Number != "")

metadata_sample_codes_fluB <- metadata_sample_codes_fluB_raw %>%
  filter(Attribute == "") %>%
  select(Blinded.Number, Tube.Label, Virus.Name, Serotype, Collection.Date, Comments, Special.Notes, Date.Sent.to.JCVI, BioProject.ID) %>%
  filter(Blinded.Number != "")

metadata_sample_codes_fluA_and_fluB <- metadata_sample_codes_fluA %>%
  bind_rows(metadata_sample_codes_fluB)

## Row bind the batch 1 and batch 2 together and then add the merge in the metadata. 

EMIT_JCVI_seq_batches_1_and_2_with_metadata <- UMD_batch_1_and_2_with_status %>%
  right_join(metadata_sample_codes_fluA_and_fluB, by = c("Blinded Number" = "Blinded.Number")) %>%
  rename(JCVI_bac_id = `JCVI Bac ID`, 
         blinded_number = `Blinded Number`,
         organism_name = `Organism Name provided by collaborator`,
         updated_organism_name = `Updated Organism Name (names that changed are in blue)`,
         JCVI_sequencing_notes = `Special Notes`,
         sequence_complete_check = `Complete/Draft`,
         sequence_complete_segment = segment,
         segment_note = note,
         UMD_sample_id = Tube.Label,
         virus_name = Virus.Name,
         subtype = Serotype,
         date_collection = Collection.Date,
         tube_location = Comments,
         sample_type = Special.Notes,
         date_sent_UMD_to_JCVI = Date.Sent.to.JCVI,
         JCVI_bioproject_id = BioProject.ID) %>%
  mutate(date_collection = parse_date_time(date_collection, orders = "d-b-y")) %>%
  mutate(date_sent_UMD_to_JCVI = parse_date_time(date_sent_UMD_to_JCVI, order = "d-b-y")) %>%
  mutate(UMD_subject_id = str_replace(UMD_sample_id, "_.*", "")) 

# Want to separate the sample_type from notes, and plan to use the ";" as the separator, but there were a few instances where there were multiple ";" in the cell, so to fix this, I plan to reprint these notes by substituting the second ";" with a hyphen.

a <- EMIT_JCVI_seq_batches_1_and_2_with_metadata[c(630, 631, 718, 719),]$sample_type
print(a)

EMIT_JCVI_seq_batches_1_and_2_with_metadata[630, ]$sample_type = "Fine aerosol (G-II Concentrated condensate); Roommate Pair Sample - Negative Control Sample - Swab tested positive, aerosol tested negative"

EMIT_JCVI_seq_batches_1_and_2_with_metadata[631, ]$sample_type = "Coarse aerosol (G-II Impactor); Roommate Pair Sample - Negative Control Sample - Swab tested postive, aerosol tested negative"

EMIT_JCVI_seq_batches_1_and_2_with_metadata[718, ]$sample_type = "Fine Aerosol (G-II Concentrated Condensate); Roommate Pair Sample - Negative Control Sample - Swab tested positive, aerosol tested negative"

EMIT_JCVI_seq_batches_1_and_2_with_metadata[719, ]$sample_type = "Fine Aerosol (G-II Concentrated Condensate); Roommate Pair Sample - Negative Control Sample - Swab tested positive, aerosol tested negative"

b <- EMIT_JCVI_seq_batches_1_and_2_with_metadata[c(630, 631, 718, 719),]$sample_type
print(b)

EMIT_JCVI_seq_batches_1_and_2_with_metadata_clean1 <- EMIT_JCVI_seq_batches_1_and_2_with_metadata %>%
  separate(sample_type, c("sample_type", "UMD_sample_notes"), sep = ";")

# clean up the sample_type names

unique(EMIT_JCVI_seq_batches_1_and_2_with_metadata_clean1$sample_type)

EMIT_JCVI_seq_batches_1_and_2_with_metadata_clean2 <- EMIT_JCVI_seq_batches_1_and_2_with_metadata_clean1 %>%
  mutate(sample_type = ifelse(sample_type == "Coarse aerosol (G-II Impactor)", "coarse", 
                              ifelse(sample_type == "Coarse Aerosol (G-II Impactor)", "coarse",
                                     ifelse(sample_type == "Fine aerosol (G-II Concentrated condensate)", "fine",
                                            ifelse(sample_type == "Fine Aerosol (G-II Concentrated Condensate)", "fine", 
                                                   ifelse(sample_type == "Nasopharyngeal Swab", "nps",
                                                          ifelse(sample_type == "Oropharyngeal Swab", "ops", sample_type))))))) 

unique(EMIT_JCVI_seq_batches_1_and_2_with_metadata_clean2$sample_type)

# clean up instances where sample_type is not clearly described

c <- EMIT_JCVI_seq_batches_1_and_2_with_metadata_clean2 %>%
  filter(sample_type == "Negative Control (NP Swab negative, aerosol positive)")
print(c)

EMIT_JCVI_seq_batches_1_and_2_with_metadata_clean3 <- EMIT_JCVI_seq_batches_1_and_2_with_metadata_clean2 %>%
  mutate(UMD_sample_notes = ifelse(sample_type == "Negative Control (NP Swab negative, aerosol positive)", "Negative Control (NP Swab negative, aerosol positive)", UMD_sample_notes)) %>%
  mutate(sample_type = ifelse(sample_type == "Negative Control (NP Swab negative, aerosol positive)", "nps", sample_type))

d <- EMIT_JCVI_seq_batches_1_and_2_with_metadata_clean3 %>%
  filter(sample_type == "Negative Control (NP Swab negative, aerosol positive)")
print(d)

write.csv(EMIT_JCVI_seq_batches_1_and_2_with_metadata_clean3, "/Users/jbueno/Box Sync/EMIT/EMIT_Data_Analysis_Jake/EMIT_UMD_Natural_Infection/Sequence_Data/JCVI/Metadata/EMIT_JCVI_seq_batches_1_and_2_with_metadata.csv")


## Checking Nancy Leung's 3-Climate Paper data ####

# It appears that Nancy may have used just the first day of G2 sampling from each participant. To test this, I will compute the number and percent of G2 positive samples, the GM and GSD, and the the range for the first day of G2 sampling for each participant, as well as for all days of G2 sampling combined (over all 3 potential study days). 

# We can filter the data to a subset where g2.run == 1 and then compute for NPS, coarse, and fine

PNAS_data_full_g2_run_1 <- PNAS_data_full %>%
  filter(g2.run == 1)

PNAS_data_full_g2_run_1_A_NPS_1_obs <- PNAS_data_full_g2_run_1 %>%
  filter(type == "A") %>%
  filter(sample.type == "Nasopharyngeal swab") %>%
  group_by(subject.id, dpo) %>%
  mutate(sample_copy_mean = mean(final.copies, na.rm = TRUE)) %>% # doing this because there were some instances where there were repeated NPS pcr assays - taking the mean of these repeats here. %>%
  distinct(subject.id, dpo, final.copies) %>%
  summarize(count = n()) %>%
  filter(count == 1) %>%
  select(-count) %>%
  left_join(PNAS_data_full_g2_run_1) %>%
  filter(type == "A") %>%
  filter(sample.type == "Nasopharyngeal swab") %>%
  select(subject.id, dpo, final.copies)
  

PNAS_data_full_g2_run_1_A_NPS_2_plus_obs <- PNAS_data_full_g2_run_1 %>%
  filter(type == "A") %>%
  filter(sample.type == "Nasopharyngeal swab") %>%
  group_by(subject.id, dpo) %>%
  mutate(sample_copy_mean = mean(final.copies, na.rm = TRUE)) %>% # doing this because there were some instances where there were repeated NPS pcr assays - taking the mean of these repeats here. %>%
  distinct(subject.id, dpo, final.copies) %>%
  summarize(count = n()) %>%
  filter(count >= 2) %>%
  select(-count) %>%
  left_join(PNAS_data_full_g2_run_1) %>%
  filter(type == "A") %>%
  filter(sample.type == "Nasopharyngeal swab") %>%
  select(subject.id, dpo, final.copies) %>%
  mutate(final.copies = if_else(is.na(final.copies), 2000*1/sqrt(2), final.copies)) %>% 
  group_by(subject.id, dpo) %>%
  mutate(final.copies = mean(final.copies)) %>%
  distinct(subject.id, dpo, final.copies)


NPS_A <- PNAS_data_full_g2_run_1_A_NPS_1_obs %>%
  bind_rows(PNAS_data_full_g2_run_1_A_NPS_2_plus_obs) %>%
  filter(!is.na(final.copies)) %>%
  mutate(ln.final.copies = log(final.copies)) %>%
  ungroup() %>%
  summarize(A_Fine_Positive_Samples_Geom_Mean = exp(mean(ln.final.copies, na.rm = TRUE)),
            A_Fine_Positive_Samples_GSD = exp(sd(ln.final.copies, na.rm = TRUE)),
            A_Fine_Positive_Samples_Max = exp(max(ln.final.copies, na.rm = TRUE)),
            A_Fine_Positive_Samples_n = n())
NPS_A$A_Fine_Positive_Samples_Geom_Mean <- 
  format(NPS_A$A_Fine_Positive_Samples_Geom_Mean, scientific = TRUE)
NPS_A$A_Fine_Positive_Samples_Max <- 
  format(NPS_A$A_Fine_Positive_Samples_Max, scientific = TRUE)

### Repeating the above for flu B using only the first day of g2 runs
PNAS_data_full_g2_run_1 <- PNAS_data_full %>%
  filter(g2.run == 1)

PNAS_data_full_g2_run_1_B_NPS_1_obs <- PNAS_data_full_g2_run_1 %>%
  filter(type == "B") %>%
  filter(sample.type == "Nasopharyngeal swab") %>%
  group_by(subject.id, dpo) %>%
  mutate(sample_copy_mean = mean(final.copies, na.rm = TRUE)) %>% # doing this because there were some instances where there were repeated NPS pcr assays - taking the mean of these repeats here. %>%
  distinct(subject.id, dpo, final.copies) %>%
  summarize(count = n()) %>%
  filter(count == 1) %>%
  select(-count) %>%
  left_join(PNAS_data_full_g2_run_1) %>%
  filter(type == "B") %>%
  filter(sample.type == "Nasopharyngeal swab") %>%
  select(subject.id, dpo, final.copies)


PNAS_data_full_g2_run_1_B_NPS_2_plus_obs <- PNAS_data_full_g2_run_1 %>%
  filter(type == "B") %>%
  filter(sample.type == "Nasopharyngeal swab") %>%
  group_by(subject.id, dpo) %>%
  mutate(sample_copy_mean = mean(final.copies, na.rm = TRUE)) %>% # doing this because there were some instances where there were repeated NPS pcr assays - taking the mean of these repeats here. %>%
  distinct(subject.id, dpo, final.copies) %>%
  summarize(count = n()) %>%
  filter(count >= 2) %>%
  select(-count) %>%
  left_join(PNAS_data_full_g2_run_1) %>%
  filter(type == "B") %>%
  filter(sample.type == "Nasopharyngeal swab") %>%
  select(subject.id, dpo, final.copies) %>%
  mutate(final.copies = if_else(is.na(final.copies), 2000*1/sqrt(2), final.copies)) %>% 
  group_by(subject.id, dpo) %>%
  mutate(final.copies = mean(final.copies)) %>%
  distinct(subject.id, dpo, final.copies)


NPS_B <- PNAS_data_full_g2_run_1_B_NPS_1_obs %>%
  bind_rows(PNAS_data_full_g2_run_1_B_NPS_2_plus_obs) %>%
  filter(!is.na(final.copies)) %>%
  mutate(ln.final.copies = log(final.copies)) %>%
  ungroup() %>%
  summarize(B_Fine_Positive_Samples_Geom_Mean = exp(mean(ln.final.copies, na.rm = TRUE)),
            B_Fine_Positive_Samples_GSD = exp(sd(ln.final.copies, na.rm = TRUE)),
            B_Fine_Positive_Samples_Max = exp(max(ln.final.copies, na.rm = TRUE)),
            B_Fine_Positive_Samples_n = n())
NPS_B$B_Fine_Positive_Samples_Geom_Mean <- 
  format(NPS_B$B_Fine_Positive_Samples_Geom_Mean, scientific = TRUE)
NPS_B$B_Fine_Positive_Samples_Max <- 
  format(NPS_B$B_Fine_Positive_Samples_Max, scientific = TRUE)





### Repeating the above but using all the data, not just the first day of g2 sampling
# For flu A
PNAS_data_full_g2_run_1_A_NPS_1_obs <- PNAS_data_full %>%
  filter(type == "A") %>%
  filter(sample.type == "Nasopharyngeal swab") %>%
  group_by(subject.id, dpo) %>%
  mutate(sample_copy_mean = mean(final.copies, na.rm = TRUE)) %>% # doing this because there were some instances where there were repeated NPS pcr assays - taking the mean of these repeats here. %>%
  distinct(subject.id, dpo, final.copies) %>%
  summarize(count = n()) %>%
  filter(count == 1) %>%
  select(-count) %>%
  left_join(PNAS_data_full) %>%
  filter(type == "A") %>%
  filter(sample.type == "Nasopharyngeal swab") %>%
  select(subject.id, dpo, final.copies)


PNAS_data_full_g2_run_1_A_NPS_2_plus_obs <- PNAS_data_full %>%
  filter(type == "A") %>%
  filter(sample.type == "Nasopharyngeal swab") %>%
  group_by(subject.id, dpo) %>%
  mutate(sample_copy_mean = mean(final.copies, na.rm = TRUE)) %>% # doing this because there were some instances where there were repeated NPS pcr assays - taking the mean of these repeats here. %>%
  distinct(subject.id, dpo, final.copies) %>%
  summarize(count = n()) %>%
  filter(count >= 2) %>%
  select(-count) %>%
  left_join(PNAS_data_full) %>%
  filter(type == "A") %>%
  filter(sample.type == "Nasopharyngeal swab") %>%
  select(subject.id, dpo, final.copies) %>%
  mutate(final.copies = if_else(is.na(final.copies), 2000*1/sqrt(2), final.copies)) %>% 
  group_by(subject.id, dpo) %>%
  mutate(final.copies = mean(final.copies)) %>%
  distinct(subject.id, dpo, final.copies)


NPS_A <- PNAS_data_full_g2_run_1_A_NPS_1_obs %>%
  bind_rows(PNAS_data_full_g2_run_1_A_NPS_2_plus_obs) %>%
  filter(!is.na(final.copies)) %>%
  mutate(ln.final.copies = log(final.copies)) %>%
  ungroup() %>%
  summarize(A_Fine_Positive_Samples_Geom_Mean = exp(mean(ln.final.copies, na.rm = TRUE)),
            A_Fine_Positive_Samples_GSD = exp(sd(ln.final.copies, na.rm = TRUE)),
            A_Fine_Positive_Samples_Max = exp(max(ln.final.copies, na.rm = TRUE)),
            A_Fine_Positive_Samples_n = n())
NPS_A$A_Fine_Positive_Samples_Geom_Mean <- 
  format(NPS_A$A_Fine_Positive_Samples_Geom_Mean, scientific = TRUE)
NPS_A$A_Fine_Positive_Samples_Max <- 
  format(NPS_A$A_Fine_Positive_Samples_Max, scientific = TRUE)

# For flu B

PNAS_data_full_g2_run_1_B_NPS_1_obs <- PNAS_data_full %>%
  filter(type == "B") %>%
  filter(sample.type == "Nasopharyngeal swab") %>%
  group_by(subject.id, dpo) %>%
  mutate(sample_copy_mean = mean(final.copies, na.rm = TRUE)) %>% # doing this because there were some instances where there were repeated NPS pcr assays - taking the mean of these repeats here. %>%
  distinct(subject.id, dpo, final.copies) %>%
  summarize(count = n()) %>%
  filter(count == 1) %>%
  select(-count) %>%
  left_join(PNAS_data_full) %>%
  filter(type == "B") %>%
  filter(sample.type == "Nasopharyngeal swab") %>%
  select(subject.id, dpo, final.copies)


PNAS_data_full_g2_run_1_B_NPS_2_plus_obs <- PNAS_data_full %>%
  filter(type == "B") %>%
  filter(sample.type == "Nasopharyngeal swab") %>%
  group_by(subject.id, dpo) %>%
  mutate(sample_copy_mean = mean(final.copies, na.rm = TRUE)) %>% # doing this because there were some instances where there were repeated NPS pcr assays - taking the mean of these repeats here. %>%
  distinct(subject.id, dpo, final.copies) %>%
  summarize(count = n()) %>%
  filter(count >= 2) %>%
  select(-count) %>%
  left_join(PNAS_data_full) %>%
  filter(type == "B") %>%
  filter(sample.type == "Nasopharyngeal swab") %>%
  select(subject.id, dpo, final.copies) %>%
  mutate(final.copies = if_else(is.na(final.copies), 2000*1/sqrt(2), final.copies)) %>% 
  group_by(subject.id, dpo) %>%
  mutate(final.copies = mean(final.copies)) %>%
  distinct(subject.id, dpo, final.copies)


NPS_B <- PNAS_data_full_g2_run_1_B_NPS_1_obs %>%
  bind_rows(PNAS_data_full_g2_run_1_B_NPS_2_plus_obs) %>%
  filter(!is.na(final.copies)) %>%
  mutate(ln.final.copies = log(final.copies)) %>%
  ungroup() %>%
  summarize(B_Fine_Positive_Samples_Geom_Mean = exp(mean(ln.final.copies, na.rm = TRUE)),
            B_Fine_Positive_Samples_GSD = exp(sd(ln.final.copies, na.rm = TRUE)),
            B_Fine_Positive_Samples_Max = exp(max(ln.final.copies, na.rm = TRUE)),
            B_Fine_Positive_Samples_n = n())
NPS_B$B_Fine_Positive_Samples_Geom_Mean <- 
  format(NPS_B$B_Fine_Positive_Samples_Geom_Mean, scientific = TRUE)
NPS_B$B_Fine_Positive_Samples_Max <- 
  format(NPS_B$B_Fine_Positive_Samples_Max, scientific = TRUE)









