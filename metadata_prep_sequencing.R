# metadata_prep_sequencing.R
# Jacob Bueno de Mesquita
# Date: July 10, 2020

# Summary

# Preparing metadata for sequencing work to be done on EMIT naturally infected (2012-2013 recruitment). 
# Sequencing was done at JCVI but in an effort to try to increase depth of coverage, samples have been transferred over the last 9 months from JCVI (Gene Tan's lab, where we had originally sent our EMIT specimens) to Mt. Sinai (Harm Van Bakel's lab).
# Harm's group requires a special metadata dictated by the CEIRS program. 
# This script automates the completion of this metadata because it would be tedious to prepare some of the variables by hand. 

# load packages ####

require(tidyverse)


# read in data ####

meta_1_1_96 <- read.csv("/Users/jbueno/Box/EMIT/EMIT_seq/5_transfer_samples_jcvi_to_mt_sinai/metadata_submission_umd_samples_from_jcvi_to_mtsinai/draft_metadata_for_ceirs_digs_submission/DPCC_Data_Standard_Template_for_Sequencing_Request_v1.1_1.csv") 

meta_2_97_192 <- read.csv("/Users/jbueno/Box/EMIT/EMIT_seq/5_transfer_samples_jcvi_to_mt_sinai/metadata_submission_umd_samples_from_jcvi_to_mtsinai/draft_metadata_for_ceirs_digs_submission/DPCC_Data_Standard_Template_for_Sequencing_Request_v1.1_2.csv") 

meta_3_193_288 <- read.csv("/Users/jbueno/Box/EMIT/EMIT_seq/5_transfer_samples_jcvi_to_mt_sinai/metadata_submission_umd_samples_from_jcvi_to_mtsinai/draft_metadata_for_ceirs_digs_submission/DPCC_Data_Standard_Template_for_Sequencing_Request_v1.1_3.csv") 

meta_4_289_325 <- read.csv("/Users/jbueno/Box/EMIT/EMIT_seq/5_transfer_samples_jcvi_to_mt_sinai/metadata_submission_umd_samples_from_jcvi_to_mtsinai/draft_metadata_for_ceirs_digs_submission/DPCC_Data_Standard_Template_for_Sequencing_Request_v1.1_4.csv") 


# add strain_name (X.3 variable)

meta_1_1_96_prep <- meta_1_1_96 %>%
  mutate(X.3 = paste("A/human/Maryland/", X, "/2012-2013", "(", X.4, ")", sep = ""))
  
meta_2_97_192_prep <- meta_2_97_192 %>%
  mutate(X.3 = paste("A/human/Maryland/", X, "/2012-2013", "(", X.4, ")", sep = ""))

meta_3_193_2886_prep <- meta_3_193_288 %>%
  mutate(X.3 = paste("A/human/Maryland/", X, "/2012-2013", "(", X.4, ")", sep = ""))

meta_4_289_325_prep <- meta_4_289_325 %>%
  mutate(X.3 = paste("A/human/Maryland/", X, "/2012-2013", "(", X.4, ")", sep = ""))

# Write out files

write.csv(meta_1_1_96_prep,
          "/Users/jbueno/Box/EMIT/EMIT_seq/5_transfer_samples_jcvi_to_mt_sinai/metadata_submission_umd_samples_from_jcvi_to_mtsinai/final_sent_metadata_for_ceirs_digs_submission/DPCC_Data_Standard_Template_for_Sequencing_Request_v1.1_1_sent.csv")

write.csv(meta_2_97_192_prep,
          "/Users/jbueno/Box/EMIT/EMIT_seq/5_transfer_samples_jcvi_to_mt_sinai/metadata_submission_umd_samples_from_jcvi_to_mtsinai/final_sent_metadata_for_ceirs_digs_submission/DPCC_Data_Standard_Template_for_Sequencing_Request_v1.1_2_sent.csv")

write.csv(meta_3_193_2886_prep,
          "/Users/jbueno/Box/EMIT/EMIT_seq/5_transfer_samples_jcvi_to_mt_sinai/metadata_submission_umd_samples_from_jcvi_to_mtsinai/final_sent_metadata_for_ceirs_digs_submission/DPCC_Data_Standard_Template_for_Sequencing_Request_v1.1_3_sent.csv")

write.csv(meta_4_289_325_prep,
          "/Users/jbueno/Box/EMIT/EMIT_seq/5_transfer_samples_jcvi_to_mt_sinai/metadata_submission_umd_samples_from_jcvi_to_mtsinai/final_sent_metadata_for_ceirs_digs_submission/DPCC_Data_Standard_Template_for_Sequencing_Request_v1.1_4_sent.csv")

