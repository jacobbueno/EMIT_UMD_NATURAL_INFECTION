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

umd_clean <- umd_emit %>%
  filter(!is.na(sample.type)) %>%
  group_by(subject.id, date.visit, sample.type) %>%
  filter(!is.na(final.copies)) %>%
  mutate(RNA_load = mean(final.copies)) %>%
  distinct(subject.id,
           sample.id,
           date.visit,
           dpo,
           sample.type, 
           RNA_load) %>%
  ungroup() %>%
  mutate(date.visit = as.Date(date.visit)) %>%
  filter(sample.type == "GII condensate NO mask" |
           sample.type == "Impactor 5 um NO mask" |
           sample.type == "Nasopharyngeal swab" |
           sample.type == "Throat swab") %>%
  mutate(sample.type = if_else(sample.type == "GII condensate NO mask", "Fine Aerosol (G-II Concentrated Condensate)", 
                               if_else(sample.type == "Impactor 5 um NO mask", "Coarse Aerosol (G-II Impactor)", if_else(sample.type == "Nasopharyngeal swab", "Nasopharyngeal Swab", "Oropharyngeal Swab"))))

# merge the RNA_load information into the metadata files

flu_a_with_rna_load <- flua %>%
  mutate(date_collection = as.Date(date_collection)) %>%
  left_join(umd_clean, by = c("study_id" = "subject.id", "date_collection" = "date.visit", "sample_type" = "sample.type")) %>%
  select(-V7) %>%
  rename(notes = V6)

flu_b_with_rna_load <- flub %>%
  mutate(date_collection = as.Date(date_collection, "%m/%d/%y")) %>%
  left_join(umd_clean, by = c("study_id" = "subject.id", "date_collection" = "date.visit", "sample_type" = "sample.type"))


## Write out the updated metadata

# write.csv(flu_a_with_rna_load, "/Users/jbueno/Box/EMIT/EMIT_seq/2_seq_emit_sent_to_jcvi_2018/updated_metadata_emit_seq_jcvi_mtsinai/umd_sent_to_jcvi_for_seq/updated_metadata_to_use/flua_emit_ceirs_seq_20200824.csv")

box_write(flu_a_with_rna_load, file_name ="flua_emit_ceirs_seq_20200923.csv", dir_id = "104645258100")



# write.csv(flu_b_with_rna_load, "/Users/jbueno/Box/EMIT/EMIT_seq/2_seq_emit_sent_to_jcvi_2018/updated_metadata_emit_seq_jcvi_mtsinai/umd_sent_to_jcvi_for_seq/updated_metadata_to_use/flub_emit_ceirs_seq_20200824.csv")

box_write(flu_b_with_rna_load, file_name ="flub_emit_ceirs_seq_20200923.csv", dir_id = "104645258100")












