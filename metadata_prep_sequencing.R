# metadata_prep_sequencing.R
# Jacob Bueno de Mesquita
# Date: July-Sept, Nov, 2020

# Summary

# Preparing metadata for sequencing work to be done on EMIT naturally infected (2012-2013 recruitment). 
# Sequencing was done at JCVI but in an effort to try to increase depth of coverage, samples have been transferred over the last 9 months from JCVI (Gene Tan's lab, where we had originally sent our EMIT specimens) to Mt. Sinai (Harm Van Bakel's lab).
# Harm's group requires a special metadata dictated by the CEIRS program. 
# This script automates the completion of this metadata because it would be tedious to prepare some of the variables by hand. 

# load packages ####

require(tidyverse)
require(boxr)

boxclient <- Sys.getenv("BOX_CLIENT_ID")
boxsecret <- Sys.getenv("BOX_CLIENT_SECRET")
box_auth()

# read in data ####

# meta_1_1_96 <- read.csv("/Users/jbueno/Box/EMIT/EMIT_seq/5_transfer_samples_jcvi_to_mt_sinai/metadata_submission_umd_samples_from_jcvi_to_mtsinai/draft_metadata_for_ceirs_digs_submission/DPCC_Data_Standard_Template_for_Sequencing_Request_v1.1_1.csv") 

meta_1_1_96 <- box_search("DPCC_Data_Standard_Template_for_Sequencing_Request_v1.1_1.csv", ancestor_folder_ids = "117680706813") %>%
  box_read_csv()


# meta_2_97_192 <- read.csv("/Users/jbueno/Box/EMIT/EMIT_seq/5_transfer_samples_jcvi_to_mt_sinai/metadata_submission_umd_samples_from_jcvi_to_mtsinai/draft_metadata_for_ceirs_digs_submission/DPCC_Data_Standard_Template_for_Sequencing_Request_v1.1_2.csv") 

meta_2_97_192 <- box_search("DPCC_Data_Standard_Template_for_Sequencing_Request_v1.1_2.csv", ancestor_folder_ids = "117680706813") %>%
  box_read_csv()


# meta_3_193_288 <- read.csv("/Users/jbueno/Box/EMIT/EMIT_seq/5_transfer_samples_jcvi_to_mt_sinai/metadata_submission_umd_samples_from_jcvi_to_mtsinai/draft_metadata_for_ceirs_digs_submission/DPCC_Data_Standard_Template_for_Sequencing_Request_v1.1_3.csv") 

meta_3_193_288 <- box_search("DPCC_Data_Standard_Template_for_Sequencing_Request_v1.1_3.csv", ancestor_folder_ids = "117680706813") %>%
  box_read_csv()


# meta_4_289_325 <- read.csv("/Users/jbueno/Box/EMIT/EMIT_seq/5_transfer_samples_jcvi_to_mt_sinai/metadata_submission_umd_samples_from_jcvi_to_mtsinai/draft_metadata_for_ceirs_digs_submission/DPCC_Data_Standard_Template_for_Sequencing_Request_v1.1_4.csv") 

meta_4_289_325 <- box_search("DPCC_Data_Standard_Template_for_Sequencing_Request_v1.1_4.csv", ancestor_folder_ids = "117680706813") %>%
  box_read_csv()

# add strain_name (X.3 variable)

meta_1_1_96_prep <- meta_1_1_96 %>%
  mutate(X.3 = paste("A/human/Maryland/", V3, "/2012", "(", V7, ")", sep = ""))
  
meta_2_97_192_prep <- meta_2_97_192 %>%
  mutate(X.3 = paste("A/human/Maryland/", V2, "/2012", "(", V6, ")", sep = ""))

meta_3_193_2886_prep <- meta_3_193_288 %>%
  mutate(X.3 = paste("A/human/Maryland/", V2, "/2012", "(", V6, ")", sep = ""))

meta_4_289_325_prep <- meta_4_289_325 %>%
  mutate(X.3 = paste("A/human/Maryland/", V2, "/2012", "(", V6, ")", sep = ""))

# Write out files

# write.csv(meta_1_1_96_prep,
#           "/Users/jbueno/Box/EMIT/EMIT_seq/5_transfer_samples_jcvi_to_mt_sinai/metadata_submission_umd_samples_from_jcvi_to_mtsinai/final_sent_metadata_for_ceirs_digs_submission/DPCC_Data_Standard_Template_for_Sequencing_Request_v1.1_1_sent.csv")

box_save(meta_1_1_96_prep, file_name ="DPCC_Data_Standard_Template_for_Sequencing_Request_v1.1_1_sent.csv", dir_id = "117681765018")



# write.csv(meta_2_97_192_prep,
#           "/Users/jbueno/Box/EMIT/EMIT_seq/5_transfer_samples_jcvi_to_mt_sinai/metadata_submission_umd_samples_from_jcvi_to_mtsinai/final_sent_metadata_for_ceirs_digs_submission/DPCC_Data_Standard_Template_for_Sequencing_Request_v1.1_2_sent.csv")

box_save(meta_2_97_192_prep, file_name ="DPCC_Data_Standard_Template_for_Sequencing_Request_v1.1_2_sent.csv", dir_id = "117681765018")



# write.csv(meta_3_193_2886_prep,
#           "/Users/jbueno/Box/EMIT/EMIT_seq/5_transfer_samples_jcvi_to_mt_sinai/metadata_submission_umd_samples_from_jcvi_to_mtsinai/final_sent_metadata_for_ceirs_digs_submission/DPCC_Data_Standard_Template_for_Sequencing_Request_v1.1_3_sent.csv")

box_save(meta_3_193_2886_prep, file_name ="DPCC_Data_Standard_Template_for_Sequencing_Request_v1.1_3_sent.csv", dir_id = "117681765018")




# write.csv(meta_4_289_325_prep,
#           "/Users/jbueno/Box/EMIT/EMIT_seq/5_transfer_samples_jcvi_to_mt_sinai/metadata_submission_umd_samples_from_jcvi_to_mtsinai/final_sent_metadata_for_ceirs_digs_submission/DPCC_Data_Standard_Template_for_Sequencing_Request_v1.1_4_sent.csv")

box_save(meta_4_289_325_prep, file_name ="DPCC_Data_Standard_Template_for_Sequencing_Request_v1.1_4_sent.csv", dir_id = "117681765018")




## Adding viral load to metadata to probe relationships between viral quantity and sequence: NOT USING TOBIT HERE BUT RATHER ELIMINATING REPLICATES BELOW LOD -- FUTURE ITERATIONS MAY PREFER TO USE TOBIT ####

# The files to which we wish to add the viral load data are: 
#"flua_emit_ceirs_seq_20200807.csv" &
#"flub_emit_ceirs_seq_20200807.csv"
# Note these files were converted to csv from xlsx format. Here we wish to iterate on the .csv versions by creating new versions of the metadata that can be used in sequence analyses. To differentiate -- we will update the date in the df title. 
# Note from sept 23, 2020: it appears that the flua dataset does not contain viral load information in it. We will have to get this information from elsewhere. Probably the flua dataset was created in response to a specific request from the CEIRS sequencing core at Mt. Sinai regarding sample names. Thus, viral load wasn't included in this dataset. It is unclear why the viral load data exists for the flub df here but not the flua df. 



# flua <- read.csv("/Users/jbueno/Box/EMIT/EMIT_seq/2_seq_emit_sent_to_jcvi_2018/updated_metadata_emit_seq_jcvi_mtsinai/umd_sent_to_jcvi_for_seq/updated_metadata_to_use/flua_emit_ceirs_seq_20200807.csv")

flua <- box_search("flua_emit_ceirs_seq_20200807.csv", ancestor_folder_ids = "104645258100") %>%
  box_read_csv()


# flub <- read.csv("/Users/jbueno/Box/EMIT/EMIT_seq/2_seq_emit_sent_to_jcvi_2018/updated_metadata_emit_seq_jcvi_mtsinai/umd_sent_to_jcvi_for_seq/updated_metadata_to_use/flub_emit_ceirs_seq_20200807.csv")

flub <- box_search("flub_emit_ceirs_seq_20200807.csv", ancestor_folder_ids = "104645258100") %>%
  box_read_csv()

# read in the "all_cases.csv" df in order to add the the viral load information to the metadata files for fluA and fluB

# umd_emit <- read.csv("/Users/jbueno/Box/EMIT/EMIT_Data_Analysis_Jake/EMIT_UMD_Natural_Infection/Curated Data/Analytical Datasets/all_cases.csv")

umd_emit <- box_search("all_cases.csv", ancestor_folder_ids = "61216893198") %>%
  box_read_csv()



# Adding on another approach to isolating the shedding data by not including RNA data for the experiment instances where there was at least one replicate below detection limit, but marking them as "discordant". See if this changes the metadata previously produced (flu_a_with_rna_load, and flu_b_with_rna_load files).

# Need to add this filter step because there are experiments testing for fluA and fluB even though often there isn't coinfection with A and B types so one of the experiments being negative is taken as negative and we examine the data for the assays that probe for the subtype where there is viral detection. The NP swab is the indicator of subtype and CDC panel was run on all of those samples. 

umd_clean <- umd_emit %>%
  filter(!is.na(sample.type)) %>%
  filter((type.inf == "H3N2" & type == "A") |
           (type.inf == "Pandemic H1" & type == "A") |
           (type.inf == "Unsubtypable A" & type == "A") |
           (type.inf == "H3N2 and PH1" & type == "A") |
           (type.inf == "H3N2 and B" & type == "A") |
           (type.inf == "H3N2 and B" & type == "B") |
           (type.inf == "B and unsubtypable A" & type == "A") |
           (type.inf == "B and unsubtypable A" & type == "B") |
           (type.inf == "B" & type == "B") ) %>%
  group_by(subject.id, date.visit, sample.type, type, Experiment) %>%
  select(subject.id,
           sample.id,
           date.visit,
           dpo,
           Experiment,
           sample.type, 
           type.inf,
           type,
           final.copies) %>%
  distinct(subject.id,
           sample.id,
           date.visit,
           dpo,
           Experiment,
           sample.type, 
           type.inf,
           type,
           final.copies) %>%
  ungroup() %>%
  mutate(date.visit = as.Date(date.visit)) %>%
  filter(sample.type == "GII condensate NO mask" |
           sample.type == "Impactor 5 um NO mask" |
           sample.type == "Nasopharyngeal swab" |
           sample.type == "Throat swab") %>%
  mutate(sample.type = if_else(sample.type == "GII condensate NO mask", "Fine", 
                               if_else(sample.type == "Impactor 5 um NO mask", "Coarse", if_else(sample.type == "Nasopharyngeal swab", "Nasopharyngeal", "Oropharyngeal"))))



umd_clean_repl_bel_lod <- umd_clean %>%
  filter(!is.na(sample.type)) %>%
  group_by(subject.id, date.visit, sample.type, Experiment) %>%
  summarize(below_lod_repl = sum(is.na(final.copies)))

umd_clean_repl_above_lod <- umd_clean %>%
  filter(!is.na(sample.type)) %>%
  group_by(subject.id, date.visit, sample.type, Experiment) %>%
  summarize(above_lod_repl = sum(!is.na(final.copies)))



# dif. between repl below and repl above lod

below_above_lod <- umd_clean_repl_bel_lod %>%
  left_join(umd_clean_repl_above_lod) %>%
  mutate(all_repl_above_lod = if_else(below_lod_repl == 0 & above_lod_repl > 0, 1, 0),
         discordant = if_else(below_lod_repl > 0 & above_lod_repl > 0, 1, 0),
         all_below_lod_or_not_run = if_else(above_lod_repl == 0, 1, 0)) %>%
  filter(!is.na(Experiment)) %>%
  left_join(umd_clean, by = c("subject.id", "date.visit", "sample.type", "Experiment")) 
# each experiment here is evaluated individually for each sample. If an experiment had all replicates tested giving a qRT-PCR value above detection limit, even if a second experiment on the same sample was discordant or showed all replicates below detection limit, the experiment with all repl above detection will be counted as having all replicates above detection and used in subsequent analysis as such.



## How many sample-experiment observations (subject-date-sample_type combinations) had discordant (qRT-PCR replicates above and below detection limit)? How many had all replicates above detection? How many experiments within a subject-date-sample_type combination had multiple experiments with the same outcome, or different outcome with respect to all replicates above, below, or discordant?

samp_exp_obs_discord <- below_above_lod %>%
  distinct(subject.id, date.visit, sample.type, discordant, .keep_all = FALSE) %>%
  filter(discordant == 1) 
nrow(samp_exp_obs_discord)

samp_exp_obs_detect_all_repl <- below_above_lod %>%
  distinct(subject.id, date.visit, sample.type, all_repl_above_lod, .keep_all = FALSE) %>%
  filter(all_repl_above_lod == 1) 
nrow(samp_exp_obs_detect_all_repl)

samp_exp_obs_detect_no_repl <- below_above_lod %>%
  distinct(subject.id, date.visit, sample.type, all_below_lod_or_not_run, .keep_all = FALSE) %>%
  filter(all_below_lod_or_not_run == 1) 
nrow(samp_exp_obs_detect_no_repl)
# This is over the group where there was at least one positive sample (subtype exists and is not negative)

# prep data
below_above_lod_clean <- below_above_lod %>%
  mutate(date.visit = as.Date(date.visit)) %>%
  filter(sample.type == "Fine" |
           sample.type == "Coarse" |
           sample.type == "Nasopharyngeal" |
           sample.type == "Oropharyngeal")



# get all replicates with qRT-PCR readings above LOD
all_repl_above_lod_df <- below_above_lod_clean %>%
  filter(all_repl_above_lod == 1) %>%
  group_by(subject.id, date.visit, sample.type, type) %>%
  mutate(RNA_load = mean(final.copies)) %>%
  ungroup() %>%
  distinct(subject.id,
           sample.id,
           date.visit,
           dpo,
           sample.type,
           type.inf,
           type,
           RNA_load) %>%
  mutate(discordant = 0)
  


# viewing the discordant replicates in full
discordant_df <- below_above_lod_clean %>%
  filter(discordant == 1) %>%
  distinct(subject.id,
           sample.id,
           date.visit,
           dpo,
           sample.type,
           type.inf,
           type,
           final.copies)  



# taking all the discordant samples (with >=1 qRT-PCR replicate above detection limit and <=1 qRT-PCR replicate below detection limit) as the mean of the replicates above detection limit and noting as discordant = 1
discordant_df <- below_above_lod_clean %>%
  filter(discordant == 1) %>%
  group_by(subject.id, date.visit, sample.type, type) %>%
  mutate(RNA_load = mean(final.copies, na.rm = TRUE)) %>%
  distinct(subject.id,
           sample.id,
           date.visit,
           dpo,
           sample.type,
           type.inf,
           type,
           RNA_load) %>% #eliminating any duplicates 
  mutate(discordant = 1)
  


# row bind the discordant samples and theose with all replicates positive
full_set <- all_repl_above_lod_df %>%
  bind_rows(discordant_df)



# # noticing issue with dpo variable, will sort that here
# 
# all_repl_above_lod_df_dpo <- all_repl_above_lod_df %>%
#   distinct(study_id, date_collection, dpo) %>%
#   mutate(dpo = as.numeric(dpo)) %>%
#   filter(!is.na(dpo)) # assuming that another sample type will provide the correct dpo for all the oropharyngeal swabs that are listed as NA, which, looks accurate since there are 60 subject_id-dates in the final seq_subset df


# flu_a_vir_load_all_rep_above_lod <- all_repl_above_lod_df %>%
#   filter(type == "A")
# 
# 
# 
# flu_b_vir_load_all_rep_above_lod <- all_repl_above_lod_df %>%
#   filter(type == "B")
# 
# 
# 
# flu_a_with_rna_load_new <- flua %>%
#   mutate(date_collection = as.Date(date_collection)) %>%
#   select(-V7) %>%
#   rename(notes = V6) %>%
#   mutate(sample_type = if_else(sample_type == "Fine Aerosol (G-II Concentrated Condensate)", "Fine", 
#                                if_else(sample_type == "Coarse Aerosol (G-II Impactor)", "Coarse", if_else(sample_type == "Nasopharyngeal Swab", "Nasopharyngeal", "Oropharyngeal Swab")))) %>%
#   left_join(flu_a_vir_load_all_rep_above_lod, by = c("study_id" = "subject.id", "date_collection" = "date.visit", "sample_type" = "sample.type"))
#   
#   
#   
#   
# flu_b_with_rna_load_new <- flub %>%
#   mutate(date_collection = as.Date(date_collection, "%m/%d/%y")) %>%
#   mutate(sample_type = if_else(sample_type == "Fine Aerosol (G-II Concentrated Condensate)", "Fine", 
#                                if_else(sample_type == "Coarse Aerosol (G-II Impactor)", "Coarse", if_else(sample_type == "Nasopharyngeal Swab", "Nasopharyngeal", "Oropharyngeal Swab")))) %>%
#   left_join(flu_b_vir_load_all_rep_above_lod, by = c("study_id" = "subject.id", "date_collection" = "date.visit", "sample_type" = "sample.type"))



# merge the RNA_load information into the metadata files

umd_clean_a <- full_set %>%
  filter(type == "A")



umd_clean_b <- full_set %>%
  filter(type == "B")



flu_a_with_rna_load <- flua %>%
  mutate(sample_type = if_else(sample_type == "Fine Aerosol (G-II Concentrated Condensate)", "Fine", 
                               if_else(sample_type == "Coarse Aerosol (G-II Impactor)", "Coarse", if_else(sample_type == "Nasopharyngeal Swab", "Nasopharyngeal", "Oropharyngeal")))) %>%
  mutate(date_collection = as.Date(date_collection)) %>%
  left_join(umd_clean_a, by = c("study_id" = "subject.id", "date_collection" = "date.visit", "sample_type" = "sample.type")) %>%
  select(-V7) %>%
  rename(notes = V6)



flu_b_with_rna_load <- flub %>%
  mutate(sample_type = if_else(sample_type == "Fine Aerosol (G-II Concentrated Condensate)", "Fine", 
                               if_else(sample_type == "Coarse Aerosol (G-II Impactor)", "Coarse", if_else(sample_type == "Nasopharyngeal Swab", "Nasopharyngeal", "Oropharyngeal")))) %>%
  mutate(date_collection = as.Date(date_collection, "%m/%d/%y")) %>%
  left_join(umd_clean_b, by = c("study_id" = "subject.id", "date_collection" = "date.visit", "sample_type" = "sample.type"))



## Write out the updated metadata

box_write(flu_a_with_rna_load, file_name ="flua_emit_ceirs_seq_20201216.csv", dir_id = "104645258100")

box_write(flu_b_with_rna_load, file_name ="flub_emit_ceirs_seq_20201216.csv", dir_id = "104645258100")



### COMAPRING WITH SEQ_SUBSET THAT WAS USED IN FINAL ANALYSIS ####

# The set that made it into the final analysis from Todd's group (e.g., had sufficient sequence data)
# Note that RNA_load variable in this df is from the earlier attempt at quantifying the RNA load in the samples (that we have since amended in this script)
seq_subset <- box_search("Merged-metadata-EMIT-samples", ancestor_folder_ids = "127774011195") %>%
  box_read() 

# Let's see how the previous version compares with this version

flu_a_b_row_bind <- flu_a_with_rna_load %>%
  bind_rows(flu_b_with_rna_load) %>%
  select(-sample.id,
         -dpo)

seq_subset_new_viral_load <- seq_subset %>%
  rename(RNA_load_old = RNA_load) %>%
  left_join(flu_a_b_row_bind, by = c("jcvi_id", "umd_sample_id", "sample_type", "study_id", "date_collection"))



# noticing issue with dpo variable, will sort that here

seq_subset_dpo <- seq_subset_new_viral_load %>%
  distinct(study_id, date_collection, dpo) %>%
  mutate(dpo = as.numeric(dpo)) %>%
  filter(!is.na(dpo)) # assuming that another sample type will provide the correct dpo for all the oropharyngeal swabs that are listed as NA, which, looks accurate since there are 60 subject_id-dates in the final seq_subset df



# noticing issue with type.inf (viral subtype) variable, will sort that here

seq_subset_inftype <- seq_subset_new_viral_load %>%
  distinct(study_id, date_collection, type.inf) %>%
  filter(!is.na(type.inf)) # assuming that the subtype is listed for a different sample type (e.g., NPS, since the problem with missing type may only affect OPS data) for the same study_id-collection-date combination



# noticing issue with type (viral infection flu A or flu B) variable, will sort that here

seq_subset_a_or_b <- seq_subset_new_viral_load %>%
  distinct(study_id, date_collection, type) %>%
  filter(!is.na(type)) # assuming that the type is listed for a different sample type (e.g., NPS, since the problem with missing type may only affect OPS data) for the same study_id-collection-date combination



# add which samples had discordant replicates

samp_exp_obs_discord_clean <- samp_exp_obs_discord %>%
  mutate(sample.type = if_else(sample.type == "GII condensate NO mask", "Fine", 
                               if_else(sample.type == "Impactor 5 um NO mask", "Coarse", if_else(sample.type == "Nasopharyngeal swab", "Nasopharyngeal", "Oropharyngeal Swab")))) %>%
  mutate(date.visit = as.Date(date.visit))



# adding the full dpo, type.inf, and type vars to the df.

seq_subset_new_viral_load_clean <- seq_subset_new_viral_load %>%
  select(-dpo,
         -type.inf,
         -type) %>%
  left_join(seq_subset_dpo) %>%
  left_join(seq_subset_inftype) %>%
  left_join(seq_subset_a_or_b) %>%
  left_join(samp_exp_obs_discord_clean, by = c("study_id" = "subject.id", "date_collection" = "date.visit", "sample_type" = "sample.type", "discordant" = "discordant")) %>%
  select(study_id,
           date_collection,
           dpo,
           sample.id,
           Sample,
           jcvi_id,
           type,
           sample_type,
           discordant,
           RNA_load_old,
           RNA_load) %>%
  distinct(study_id,
           date_collection,
           dpo,
           sample.id,
           Sample,
           jcvi_id,
           type,
           sample_type,
           discordant,
           RNA_load_old,
           RNA_load) #eliminating duplicates

