# sequence_analysis.R
# Jacob Bueno de Mesquita
# Date: December 2020

# Summary

# Working with Todd Treangen's team on the sequence analysis that they are preparing using the EMIT sequences from NP, OP, and aerosols, initially sequenced at JCVI, with a subset resequenced at Mt. Sinai.
# This script was initiated with the purpose of creating an table describing population composition socio demographics/health history/symptoms during acute infection. 
# Since the samples in this sequencing analysis represent a subset of the samples analyzed and reported on in Yan et al., 2018 (PNAS) we will provide a new version of the Table 1 from that paper.
# May be interesting to look to see if there is any systematic difference in any of the variables (particularly viral shedding quantity) between the full PNAS dataset and the data for which samples were available for sequencing and could be found in the chest freezer, and resulted in producing acceptable sequence data to move forward with analysis. 

# load packages ####

require(tidyverse)
require(boxr)
require(kableExtra)
require(arsenal)
require(gridExtra)

boxclient <- Sys.getenv("BOX_CLIENT_ID")
boxsecret <- Sys.getenv("BOX_CLIENT_SECRET")
box_auth()

set.seed(1202)


# load data ####



# load the set of data that contains at a minimum, subject id, sample type, collection date, and day post symptom onset that Todd's team included in the analysis

seq_subset <- box_search("Merged-metadata-EMIT-samples", ancestor_folder_ids = "127774011195") %>%
  box_read() 



# load the full EMIT dataset that contains the symptom data of interest

emit_nat <- box_search("PNAS_data_full.csv", ancestor_folder_ids = "61216893198") %>% 
  box_read()



# prepare to merge these data files ####
# eliminate the viral load portion of the emit data -- viral load is already in the seq_subset df

emit_nat_clean <- emit_nat %>% 
  mutate(up_resp_sc = (nose_run + nose_stuf + sneeze + throat_sr + earache),
         low_resp_sc = (chest_tight + sob + cough),
         systemic_sc = (malaise + headache + mj_ache + sw_fever_chill)) %>%
  mutate(fluvac_cur = as.numeric(fluvac_cur)) %>%
  mutate(vax_bothyear = as.numeric(vax_bothyear)) %>%
  mutate(fluvac_cur = if_else(is.na(fluvac_cur), 0, fluvac_cur)) %>%
  mutate(vax_bothyear = if_else(is.na(vax_bothyear), 0, vax_bothyear)) %>%
  mutate(cens = if_else(is.na(final.copies), 1, 0)) %>%
  mutate(final.copies.lod = if_else(is.na(final.copies), 500, final.copies)) %>%
  arrange(subject.id, 
          date.visit) %>%
  mutate(log10_final_copies_lod = log10(final.copies.lod)) %>%
  mutate(date.visit = as.Date(date.visit)) %>%
  mutate(total_sx_score = up_resp_sc + low_resp_sc + systemic_sc) %>%
  distinct(subject.id,
           date.visit,
           cough_number,
           sneeze_number,
           type.inf,
           body_temp,
           nose_run,
           nose_stuf,
           sneeze,
           throat_sr,
           earache,
           malaise,
           headache,
           mj_ache,
           sw_fever_chill,
           lymph_node,
           chest_tight,
           sob,
           cough,
           asthma,
           Smoker,
           sex,
           age,
           BMI,
           fluvac_cur,
           vax_bothyear,
           up_resp_sc,
           low_resp_sc,
           systemic_sc,
           total_sx_score,
           anitviral_24h,
           .keep_all = FALSE
  )


# merge seq subset with sx and sociodem/health hist. data ####

seq_subset_tab1 <- seq_subset %>%
  left_join(emit_nat_clean, by = c("study_id" = "subject.id", "date_collection" = "date.visit"))
 

# preparing table for the total participants with sequence data on any of the 4 sample types...
# vars that remain constant regarding of sample day (age, sex, bmi, vax history, health history)
# later we will get to symptom data that does change between days over the course of sampling during illness


# determine how many subjects ever had febrile illness >37.9C

ever_fever <- seq_subset_tab1 %>%
  mutate(fever_37_9 = if_else(body_temp > 37.9, 1, 0)) %>%
  filter(!is.na(fever_37_9)) %>%
  group_by(study_id) %>%
  summarise(count_fever = sum(fever_37_9)) %>%
  mutate(ever_fever = if_else(count_fever >= 1, 1, 0)) %>%
  select(-count_fever) %>%
  ungroup()


# determine how many subjects were asymptomatic on all days of sx observation (range 1-3 days)

asympt <- seq_subset_tab1 %>%
  mutate(total_sx_score = up_resp_sc + low_resp_sc + systemic_sc) %>%
  mutate(sympt_instances = if_else(total_sx_score > 0, 1, 0)) %>%
  filter(!is.na(sympt_instances)) %>%
  group_by(study_id) %>%
  summarise(sympt_inst_count = sum(sympt_instances)) %>%
  mutate(asymptomatic = if_else(sympt_inst_count > 0, 0, 1)) %>%
  select(-sympt_inst_count) %>%
  ungroup()



# determine how many subjects have data from 1, 2, 3, or 4 types of samples

samp_type_n <- seq_subset_tab1 %>%
  distinct(study_id,
           sample_type, 
           .keep_all = FALSE) %>%
  group_by(sample_type) %>%
  summarise(n = n())
  
samp_type_sid <- seq_subset_tab1 %>%
  distinct(study_id,
           sample_type, 
           .keep_all = FALSE) %>%
  group_by(study_id) %>%
  summarise(n = n()) %>%
  group_by(n) %>%
  summarise(n_samp_types = n())

n_samp_types_per_sid <- seq_subset_tab1 %>%
  distinct(study_id,
           sample_type, 
           .keep_all = FALSE) %>%
  group_by(study_id) %>%
  summarise(n_sample_types = n())



# determine how many subjects have data from multiple days of sampling

samp_type_n <- seq_subset_tab1 %>%
  distinct(study_id,
           date_collection, 
           .keep_all = FALSE) %>%
  group_by(study_id) %>%
  summarise(n_samp_inst = n())

summary_samp_type_n <- seq_subset_tab1 %>%
  distinct(study_id,
           date_collection, 
           .keep_all = FALSE) %>%
  group_by(study_id) %>%
  summarise(n = n()) %>%
  group_by(n) %>%
  summarise(summary_n_samp_instances = n())



# now prepping data for table creation

seq_sid <- seq_subset_tab1 %>%
  distinct(study_id, .keep_all = T) %>%
  left_join(ever_fever) %>%
  left_join(n_samp_types_per_sid) %>%
  left_join(samp_type_n) %>%
  left_join(asympt) %>%
  mutate(sex = as.factor(sex),
         asthma = as.factor(asthma),
         Smoker = as.factor(Smoker),
         fluvac_cur = as.factor(fluvac_cur),
         vax_bothyear = as.factor(vax_bothyear),
         ever_fever = as.factor(ever_fever),
         n_sample_types = as.factor(n_sample_types),
         n_samp_inst = as.factor(n_samp_inst),
         asymptomatic = as.factor(asymptomatic))



# demographics/single obs per subject table

dem_summary <- tableby(  ~age +
                         sex +
                         BMI +
                         asthma +
                         Smoker +
                         fluvac_cur +
                         vax_bothyear +
                           ever_fever +
                           asymptomatic +
                           n_sample_types +
                           n_samp_inst,
                       data = seq_sid,
                       control = tableby.control(numeric.stats = "medianq1q3", digits = 1L, total = FALSE, test = TRUE),
                       #cat.simplify = TRUE,
                       numeric.simplify = TRUE)
dem_sum <- summary(dem_summary, text = TRUE)

kable(dem_sum, align = "c") %>%
  kable_styling(c("striped", "bordered", "hover", "condensed", "responsive"), 
                full_width = FALSE, position = "center", font_size = 10) %>%
  kable_styling(full_width = F) %>%
  column_spec(1, bold = T)



## Now getting to the symptom variables table prepared


seq_subset_tab1_prep <- seq_subset_tab1 %>%
  mutate(nose_run = as.numeric(nose_run),
         nose_stuf = as.numeric(nose_stuf),
         sneeze = as.numeric(sneeze),
         throat_sr = as.numeric(throat_sr),
         earache = as.numeric(earache),
         malaise = as.numeric(malaise),
         headache = as.numeric(headache),
         mj_ache = as.numeric(mj_ache),
         sw_fever_chill = as.numeric(sw_fever_chill),
         chest_tight = as.numeric(chest_tight),
         sob = as.numeric(sob),
         cough = as.numeric(cough),
         antiviral_24h = as.factor(antiviral_24h))


sx_summary <- tableby( ~ 
                        body_temp +
                         cough_number +
                         sneeze_number +
                         nose_run +
                         nose_stuf +
                         sneeze +
                         throat_sr +
                         earache +
                         malaise +
                         headache +
                         mj_ache +
                         sw_fever_chill +
                         chest_tight +
                         sob +
                         cough +
                        up_resp_sc + 
                        low_resp_sc + 
                        systemic_sc +  
                        total_sx_score +
                         anitviral_24h,
                      data = seq_subset_tab1_prep,
                      control = tableby.control(numeric.stats = "medianq1q3", digits = 1L, total = FALSE, test = TRUE, numeric.test="kwt"),
                      #cat.simplify = TRUE,
                      numeric.simplify = TRUE
                      )
summary(sx_summary, text = TRUE)

sx_sum <- summary(sx_summary, text = TRUE)

kable(sx_sum, align = "c") %>%
  kable_styling(c("striped", "bordered", "hover", "condensed", "responsive"), 
                full_width = FALSE, position = "center", font_size = 10) %>%
  kable_styling(full_width = F) %>%
  column_spec(1, bold = T)








