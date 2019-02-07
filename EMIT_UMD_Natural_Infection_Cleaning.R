# EMIT UMD Natural Infection Study Data Curation - cleaning raw data to produce cleaned spreadsheets
# Program Objective: Take the datasets identified as critical, clean them, and later merge to form curated one or more curated datasets
# Author: Jacob Bueno de Mesquita using material from Jing Yan and Don Milton
# Date: December 14, 2018 - February 2019
# Summary:

#### Load required packages and set working directory ####

library(tidyverse)
library(RcppRoll)
library(readxl)
library(knitr)
library(data.table)
library(lubridate)
library(arsenal)
library(lme4)

setwd("/Users/jbueno/Box Sync/EMIT/EMIT_Data_Analysis_Jake/EMIT_UMD_Natural_Infection")

sessionInfo() # for reproducibility

# Now pasting code from Jing Yan and Don Milton that was used in previous work on the EMIT UMD data and working with it to make sure that it functions like it was intended and to reproduce the final dataset that was used in the PNAS analysis.
# The goal here is to review their script and improve the clarity

#### **** Using Script: Jing Yan and Dr. Milton's "Merge 1-3-update.R" **** ####

##
# Original file information:

# From Jing Yan & Don Milton; January 20-25, 2016
# Purpose: follow the data analysis plan in folder EMIT_Data_Analysis (described below)

# Purpose, input and output files: 
#   1.	Read in Clinical (encounter and questionnaire) data;
#       a.	Count rows and columns check that there are no missing required fields (e.g. date of visit, subject_id). 
#           Print number of rows total and numbers for eachredcap_event_name. 
#       b.	Create a dataframe with visits (1,2,3) and another df with g2_run_1,2,3; then if date_visit=”” and redcap_event_name = “g2_run_1” 
#         then date_visit = date_g2_1; next create an indicator for visit(1,2,3) called visit_num with values (1,2, or 3) and an 
#         indicator called g2_run with values (1,2, or3) in the respective df; then merge these df on date_visit and subj_id. 
#         The result of this merge is that we get one record per subj_id and date. Create an enrolled indicator: if g2_run is.na 
#         then enrolled=FALSE, if g2_run = (1, 2, or 3) then enrolled = TRUE. Print out numbers of enrolled for each of 
#         visit_num = 1, 2 and 3 (i.e. enrolled 1st, 2nd, and 3rd screening visits).
#     c.	Keep just: subj_id, date_visit, visit_num, g2_num, enrolled indicator -> sum_clinical with variable names:
#         "field_subj_id" "date_visit"    "visit_name"    "visit_num"     "g2_name"       "g2_run"        "enrolled"      "clinical.i" 
# 2.	Read in G-II_Log:
#     a.	Count rows and columns check that there are no missing required fields (e.g. date of visit, subject_id). 
#     b.	Create indicator g2_coll_num = 1, 2 or 3 based on redcap_event_name and rename start_dt to date_visit.
#     c.	Print out numbers of rows for each redcap_event_name 
#     d.  Convert date_visit from char to date format & create indicator for g2 log record (g2lm.i)
#     d.	Keep: subj_id, date_visit, g2_coll_num, and g2lm.i -> g2_log_min with variable names:
#         "field_subj_id" "date_visit"    "g2_coll_num"   "g2lm.i""field_subj_id" "date_visit"    "g2_coll_num"   "g2lm.i" 
# 3.	Merge 1:  merge sum_clinical and g2_log_min
#     a.	Merge by subj_id and date_visit
#     b.	Print data checks (e.g. number of visits, subjects, enrolled, etc.)
#     c.	Check that all records that are marked as “enrolled” have g2_coll_num that is not na otherwise delete that extra visit 
#         (e.g. 69). Print table of number of obs by number of visits, etc. as data checks.
#     d.	Output dataframe with one obs per subj_id and date -> merge1 with variable names: 
#         "field_subj_id" "date_visit"    "g2_coll_num"   "enrolled"      "visit_num"     "g2_run"        "clinical.i"    "g2lm.i"        "merge1.i"
# 4.	Read Field Db 
#     a.	Check for empty row etc. 
#     b.	Delete empty rows, create indicator for record present in Field Db (field.db1.i)
#     c.	Print number of obs and number of obs by sample type
#     d.  Convert date_visit from char to date data type
#     d.  Output file -> field.db1 with variable names:
#         "field_subj_id" "sample_id"     "date_visit"      "sample_type"   "field.db1.i"
# 5.  Merge 2: merge merge1 and field.db1 (by field_subj_id and date_visit)
#     a.  Use the enrollment indicator to identify & remove G-II samples that were not collected but included in Field Db. 
#     b.  Data checks for numbers of rows, IDs, etc.
#     c.  Output dataframe with one obs per sample_id -> merge2
# 6.	Read UMD Samples 2013 from Redcap:
#     a.	Check for empty rows etc. delete empty rows with no date_visit or subj_id
#     b.	Print number of obs and number of obs by sample type.
#     c.	Separate collection_1, assay_1, and assay_2 
#     d.  Pull out the passage and focus assays into separate dataframes (drop pcr variables)
#     d.	Compute passage and focus assay results
#     e.	Merge by sample_id to get one row per sample with culture results. 
#     f.	Keep focus and passage variables, date_visit, sample type, subj_id, sample_id, and create indicator variables.
#     g.	Print number of obs and number of obs by sample type
#     h.	Do data checks for problems (e.g. impactors that were cultured or culture results for samples without a sample type)
#     i.	List problem samples
#     j.  Output dataframe with one obs per sample -> culture_results
# 7.	Merge 3: Use merge2 and culture_results
#     a.	By sample_id and date_visit.
#     b.	Check that merge worked by count of rows and columns and obs by sample type and check for empty variables that 
#         should have values. Check that sample types match. Print data checks
#     c.	Drop obs from field db that do not have corresponding values in the Redcap sample log. 
#     d.	Output dataframe with naone obs per sample_id -> samples.cc

##

#### READ in and work with CLINICAL DATABASE ####

clinical_in_file <- 'UMD_Raw_Data/REDCAP/EMITClinicalUMD2013.csv'
clinical_umd <- read.csv(clinical_in_file)

# Clean up the raw data just a little 
# clinical_umd$date_indx_on <- as.Date(clinical_umd$date_indx_on, format = "%m/%d/%y")
# clinical_umd$date_indx_visit <- as.Date(clinical_umd$date_indx_visit, format = "%m/%d/%y")
# clinical_umd$date_ref <- as.Date(clinical_umd$date_ref, format = "%m/%d/%y")
# clinical_umd$date_enroll <- as.Date(clinical_umd$date_enroll, format = "%m/%d/%y")
# clinical_umd$date_on_sx <- as.Date(clinical_umd$date_on_sx, format = "%m/%d/%y")
# clinical_umd$date_visit <- as.Date(clinical_umd$date_visit, format = "%m/%d/%y")
# clinical_umd$date_fever_on <- as.Date(clinical_umd$date_fever_on, format = "%m/%d/%y")
# clinical_umd$date_cough <- as.Date(clinical_umd$date_cough, format = "%m/%d/%y")
# clinical_umd$date_g2_1 <- as.Date(clinical_umd$date_g2_1, format = "%m/%d/%y")

# Let's produce some summary information about this clinical_umd df

print(nrow(clinical_umd))
print(ncol(clinical_umd))

print(sum(clinical_umd$redcap_event_name == 'visit_1_part_a_arm_1'))
print(sum(clinical_umd$redcap_event_name == 'screen_visit_2_arm_1'))
print(sum(clinical_umd$redcap_event_name == 'screen_visit_3_arm_1'))
print(sum(clinical_umd$redcap_event_name == 'visit_1_part_a_arm_1') +
        sum(clinical_umd$redcap_event_name == 'screen_visit_2_arm_1') +
        sum(clinical_umd$redcap_event_name == 'screen_visit_3_arm_1'))

print(sum(clinical_umd$redcap_event_name == 'g2_run_1_arm_1'))
print(sum(clinical_umd$redcap_event_name == 'g2_run_2_arm_1'))
print(sum(clinical_umd$redcap_event_name == 'g2_run_3_arm_1'))
print(sum(clinical_umd$redcap_event_name == 'g2_run_1_arm_1') +
        sum(clinical_umd$redcap_event_name == 'g2_run_2_arm_1') +
        sum(clinical_umd$redcap_event_name == 'g2_run_3_arm_1'))

print(addmargins(with(clinical_umd, table(redcap_event_name, exclude = c()))))

# Note that one subject was enrolled twice!
print(select(filter(clinical_umd, field_subj_id == 47 | field_subj_id == 187),
             field_subj_id, date_visit, redcap_event_name, date_g2_1, rapid_flu___1, rapid_flu___2, spec_note))
# This means that there was actually one less unique invidivual than we have unique subjects ids.
# We note the above and treat subject IDs as person-illness-episodes, not persons.

# Clinical data split into screening visits and g2 runs and remerged to get one row per encounter date
clinical_min <- clinical_umd %>% 
  select(field_subj_id, redcap_event_name, date_visit, date_g2_1, rapid_flu___3, rapid_flu_loc, body_temp, date_on_sx)

clinical_visit <- clinical_min %>% 
  filter(grepl('visit', redcap_event_name))

clinical_g2 <- clinical_min %>% 
  filter(grepl('^g2_run', redcap_event_name))

clinical_visit$visit_num <- ifelse(clinical_visit$redcap_event_name == 'visit_1_part_a_arm_1', 1, 
                                   ifelse(clinical_visit$redcap_event_name == 'screen_visit_2_arm_1', 2, 3)) 

clinical_g2$g2_run <- ifelse(clinical_g2$redcap_event_name == 'g2_run_1_arm_1', 1, 
                             ifelse(clinical_g2$redcap_event_name == 'g2_run_2_arm_1', 2, 3))

# Find the G2 first time sample collection data
clinical_g2_1 <- clinical_g2 %>% 
  filter(grepl('^g2_run_1', redcap_event_name)) %>%
  select(field_subj_id, redcap_event_name, date_g2_1, g2_run) %>%
  rename(date_visit = date_g2_1)

# Find the G2 2nd and 3rd time sample collection data
clinical_g2_23 <- clinical_g2 %>%
  filter(!grepl('^g2_run_1', redcap_event_name)) %>% 
  select(field_subj_id, redcap_event_name, date_visit, g2_run)

# Merge the G2 sample collection visits (1, 2, and 3) together into the clinical_g2 df
clinical_g2 <- merge(clinical_g2_1, clinical_g2_23, 
                     c('field_subj_id', 'date_visit', 'redcap_event_name', 'g2_run'), all = TRUE)

sum_clinical <- merge(select(clinical_visit, -contains("date_g2_1")), clinical_g2, 
                      c('field_subj_id', 'date_visit'), all = TRUE)

sum_clinical$enrolled <- ifelse(!is.na(sum_clinical$g2_run), TRUE, FALSE)

sum_clinical <- sum_clinical %>% 
  rename(visit_name = redcap_event_name.x, g2_name = redcap_event_name.y)

sum_clinical$g2_run <- with(sum_clinical, ifelse(is.na(g2_run), 0, g2_run)) 

sum_clinical$visit_num <- with(sum_clinical, ifelse(is.na(visit_num), 999, visit_num))

sum_clinical$clinical.i <- TRUE #indicator for presence of record in sum_clinical

# Number of rows in summary data (total number of unique encouters)
nrow(sum_clinical)

# Tabulations for data checks
print(addmargins(with(sum_clinical, table(visit_name, visit_num, exclude = c()))))
print(addmargins(with(sum_clinical, table(g2_name, g2_run, exclude = c()))))
print(addmargins(with(sum_clinical, table(enrolled, clinical.i, exclude = c()))))
print(addmargins(with(sum_clinical, table(visit_num, g2_run, exclude = c()))))

# Total number of g2 runs (sum of runs 1, 2, and 3) according to initial clinical data
sum(sum_clinical$g2_run > 0)

# Fix the format of the date variable in sum_clinical
sum_clinical$date_visit <- as.Date(as.character(sum_clinical$date_visit), format = "%m/%d/%y")
sum_clinical$date_on_sx <- as.Date(as.character(sum_clinical$date_on_sx), format = "%m/%d/%y")
sum_clinical$date_on_sx <- as.factor(sum_clinical$date_on_sx)
print(head(tbl_df(sum_clinical)))

sum_clinical_subjectID_count <- sum_clinical %>%
  group_by(field_subj_id) %>%
  summarise(count = n())
# This shows that this df has information on all 355 screened participants

#### READ in and work with G2 LOG DATA ####

g2_in_file <- 'UMD_Raw_Data/GII/EMITGIILogUMD2013.csv'
g2_log <- read.csv(g2_in_file)

# Check the number of sampling instances here
g2_log_sampling_instance_check <- g2_log %>%
  distinct(subject_id, start_dt)
# This yields 277 but this is before the data clean. 

print(nrow(g2_log))
print(ncol(g2_log))

# Date Entry Error Correction
# Subject_id 284 g2 collection_2_arm_1 was entered as 2013-03-17 but baseline was on 2013-02-16 and collection_3 was 2013-02-18.
# Therefore recode collection_2 date to February from March (i.e. to 2013-02-17).
print(select(filter(g2_log, subject_id == 284), subject_id, redcap_event_name, start_dt))
g2_log$start_dt[which(g2_log$subject_id == 284 & g2_log$start_dt == '2013-03-17', arr.ind = TRUE)] <- '2013-02-17'

# Number of subject(s) without a start_dt
print(sum((g2_log$start_dt) == ""))

id <- as.integer(select(filter(g2_log, start_dt == ""), subject_id)) 

# G-II log for cases without start date
print(tbl_df(select(filter(g2_log, subject_id==id), 
                    subject_id, 
                    redcap_event_name, 
                    start_dt, 
                    g2_unit,operator,
                    chiller_t1,
                    subj_min)))

# Clinical data for cases without start date 
print(tbl_df(filter(sum_clinical, field_subj_id == id)))

# Data Entry Error Correction
# Remove empty record identified as having no start_dt, subject_id:
g2_log$subject_id[which(g2_log$start_dt == "", arr.ind = TRUE)]
g2_log <- g2_log %>%
  filter(!(g2_log$start_dt) == "")

# Number of rows in g2_log data
print(nrow(g2_log))
# Number of cols in g2_log data
print(ncol(g2_log))

# Number of subjects with a coll_arm_1 
print(sum(g2_log$redcap_event_name == 'baseline_and_colle_arm_1'))

# Number of subjects with a coll_2 arm
print(sum(g2_log$redcap_event_name == 'collection_2_arm_1'))

# Number of subjects with a coll_3 arm
print(sum(g2_log$redcap_event_name == 'collection_3_arm_1'))

g2_log <- g2_log %>%
  rename(date_visit = start_dt)

g2_log_min <- g2_log %>% 
  select(subject_id, redcap_event_name, date_visit, subj_min, cough_number, sneeze_number)
# Not sure if we need to keep the cough_count and subject_min data in here.
# However, note that at this point much of the cough data is missing 
# These mising cough values can be filled in based on the audio recordings.
# Adding in the audio recording cough data is done elsewhere, but need to find where?
# Here we kept in the cough_number variable but it appears to be missing some of the values that have audio and haven't yet been added into the df.

g2_log_min$g2_coll_num <- ifelse(g2_log_min$redcap_event_name == 'baseline_and_colle_arm_1', 1, 
                              ifelse(g2_log_min$redcap_event_name == 'collection_2_arm_1', 2, 3)) 

# Numbers of subjects by g2 collection event
print(ftable(addmargins(with(g2_log_min, table(redcap_event_name, g2_coll_num, exclude = c())))))
# This is an important print out of the number of G2 collection events by visit number.

g2_log_min <- g2_log_min %>% 
  select(subject_id, date_visit, g2_coll_num, subj_min, cough_number, sneeze_number) %>%
  rename(field_subj_id = subject_id)

g2_log_min$g2lm.i <- TRUE #indicator for preseence of record in g2_log_min

g2_log_min$date_visit <- as.Date(g2_log_min$date_visit)

# Check the variable names in g2_log_min
print(head(tbl_df(g2_log_min)))

g2_log_min_subjectID_count <- g2_log_min %>%
  group_by(field_subj_id) %>%
  summarise (count = n())
# 178 subjects

# Check the number of sampling instances here
g2_log_sampling_instance_check <- g2_log_min %>%
  distinct(field_subj_id, date_visit)
print(nrow(g2_log_sampling_instance_check))
# This shows that this df has information on the 276 instances of Gii runs for the 178 enrolled participants

#### MERGE CLINICAL AND G2-LOG DATA ####

# Here we will merge the sum_clinical df and the g2_log_min df
# These dfs were manipulated in the previous two sections of code in this script ...
# ... in preparation for this merging step.

merge1 <- merge(sum_clinical, g2_log_min, by = c('field_subj_id', 'date_visit'), all=TRUE)

# Merge 1 check the dimensions agains the source dfs
print(ftable(addmargins(with(merge1, table(clinical.i, g2lm.i, exclude = c())))))

# Select the variables of importance and their order. 
merge1 <- merge1 %>% 
  select(field_subj_id, 
         date_visit, 
         date_on_sx,
         g2_coll_num, 
         enrolled, 
         visit_num, 
         g2_run, 
         clinical.i, 
         g2lm.i, 
         rapid_flu___3, 
         rapid_flu_loc, 
         body_temp, 
         cough_number,
         sneeze_number)

# Do the g2_coll_num match with the g2_run? They should match.
identical(merge1$g2_coll_num, merge1$g2_run)
# However they don't match! The next lines will address the discrepancies

# If there is no g2 run then the collection number (g2_coll_num) and run number (g2_run) ...
# ... are set to zero here to match the g2_run variable.
merge1$g2_coll_num[is.na(merge1$g2_coll_num)] <- 0
merge1$g2_run[is.na(merge1$g2_run)] <- 0
merge1$indicator <- ifelse(merge1$g2_coll_num == merge1$g2_run, 1, 0) 
i <- which(merge1$indicator == 0, arr.ind = TRUE)

# Since g2_coll_num and g2_run don't match, find the unmatched subject, ...
# ... and show which row and column the subject(s) located:
i
# Subject that does not match field_subj_id
merge1$field_subj_id[i]
# Initiate data editing to resolve this discrepancy.

## Data Editing ##

# Subject 69 has data that doesn't match in the g2_coll_num and g2_run variables (which should be the same)
# The g2_coll_num variable comes from the clinical database df and the g2_run var comes from the g2_log df. 

merge1$field_subj_id[which(merge1$indicator == 0, arr.ind = TRUE)]

# Subject 69 had a second g2 visit in clinical data, but did not provide sample - 
# See REDCap comments for details - removed from final analysis data sets
merge1 <- merge1 %>% 
  filter(indicator == 1) %>% 
  select(-indicator)

# Recheck after delection: Merge 1 record sources
print(ftable(addmargins(with(merge1, table(clinical.i, g2lm.i, exclude = c())))))

# Number of rows in edited merged data
print(nrow(merge1))

# Number of subjects with a first visit
print(sum(merge1$visit_num == 1, na.rm = TRUE))

# Number of subjects with a second screen
print(sum(merge1$visit_num == 2, na.rm = TRUE))

# Number of subjects with a third screen
print(sum(merge1$visit_num == 3, na.rm = TRUE))

# Number of unique subjects
print(length(unique(merge1$field_subj_id)))

# Number of screening visits
sum(merge1$visit_num == 1, na.rm = TRUE) + 
  sum(merge1$visit_num == 2, na.rm = TRUE) + 
  sum(merge1$visit_num == 3, na.rm = TRUE)

# Number of subjects with a 1st g2 run ...
# g2_num:
print(sum(merge1$g2_run == 1))
# g2_coll_num:
print(sum(merge1$g2_coll_num == 1))

# Number of subjects with a 2nd g2 run ...
# g2_num:
print(sum(merge1$g2_run == 2))
# g2_coll_num:
print(sum(merge1$g2_coll_num == 2))

# Number of subjects with a 3rd g2 run
# g2_num:
print(sum(merge1$g2_run == 3))
# g2_coll_num: 
print(sum(merge1$g2_coll_num == 3))

# Total number of g2 runs
# based on g2_num:
print(sum(!merge1$g2_run == 0)) 
# based on g2_coll_num:
print(sum(!merge1$g2_coll_num == 0))

# Total number of screenings without a g2 run whether or not later enrolled ...
# ... based on g2_num
print(sum(merge1$g2_run == 0))
# ... based on g2_coll_num
print(sum(merge1$g2_coll_num == 0))

t1 <- sum(!merge1$g2_run == 0) + sum(merge1$g2_run == 0)
t2 <- sum(!merge1$g2_coll_num == 0) + sum(merge1$g2_coll_num == 0)

# Total number of encounters ... 
# ... based on g2_run:
t1
# ... based on g2_coll_num
t2

merge1$merge1.i <- T #indicator for record in merge1

# Cross tab tables of data in merge1
print(ftable(addmargins(with(merge1, table(visit_num, g2_run, g2_coll_num, exclude = c())))))

# Variables in merge1
print(head(tbl_df(merge1)))

#### READ in and work with the FIELD SAMPLE DATABASE ####

# Input Field Sample Data
field_db_in_file <- 'UMD_Raw_Data/EMIT UMD Field_db/field_db.csv'
field.db <- read.csv(field_db_in_file, as.is = T)

print(head(tbl_df(field.db)))

field.db1 <- field.db %>% 
  select(SUBJECT_IDENTIFIER, SAMPLE_ID, COLLECTION_DT, TYPE_NAME) %>%
  rename(field_subj_id = SUBJECT_IDENTIFIER) %>%
  rename(sample_id = SAMPLE_ID) %>%
  rename(date_visit = COLLECTION_DT) %>%
  rename(sample_type = TYPE_NAME)
field.db1$field.db1.i <- TRUE #indicator that data is in field.db1

# Input field sample field.db file
field_db_in_file

# Number of rows in field sample database
print(nrow(field.db1))
# Number of columns in database (selected columns)
ncol(field.db1)

# Tabluation of number of rows by sample type in field database
print(addmargins(with(field.db1, table(sample_type, field.db1.i, exclude = c()))))

# Number of rows that have a missing sample_type
print(sum(field.db1$sample_type == ""))

# Number of rows that have a missing date_visit
print(sum(field.db1$date_visit == ""))

# Number of rows that have a NA for date_visit
print(sum(is.na(field.db1$date_visit)))

# Number of subjects that have a missing sample_id
print(sum(field.db1$sample_id == ""))

# Number of subjects that have a NA for sample_id
print(sum(is.na(field.db1$sample_id)))

field.db1$date_visit <- as.Date(field.db1$date_visit, format = "%m/%d/%Y")

## Data Editing ##

# Field_subj_id 225 was moved to field_subj_id 250 in clinical database because second screening visit was ...
# ... erroneously given a new ID number. However the samples are still shown in the field database as 225.
# Therefore, I am recoding the subject id to 250 but leaving the sample_id as 225_x -- at least for now

field.db1$field_subj_id <- with(field.db1, ifelse(field_subj_id == 225, 250, field_subj_id))
print(tbl_df(filter(field.db1, field_subj_id == 250)))

# Recode date_visit subj 10 sample_id = 10_6 to be correct date of 2012-12-05

# orginal data
print(tbl_df(filter(field.db1, field_subj_id == 10)))

field.db1$date_visit <- with(field.db1, 
                             ifelse(field_subj_id == 10 & date_visit == as.Date("2012-12-07"), 
                                    as.Date("2012-12-05"), date_visit))
field.db1$date_visit <- as.Date(field.db1$date_visit, format = "%Y-%m-%d", origin = "1970-01-01")

# recoded data
print(tbl_df(filter(field.db1, field_subj_id == 10)))

# Recode date_visit subj 12 sample_id = 12_6 to the correct date of 2012-12-10

# orginal data
print(tbl_df(filter(field.db1, field_subj_id == 12)))
 
field.db1$date_visit <- with(field.db1, 
                             ifelse(field_subj_id == 12 & date_visit == as.Date("2013-02-08"), as.Date("2012-12-10"), date_visit))
field.db1$date_visit <- as.Date(field.db1$date_visit, format = "%Y-%m-%d", origin = "1970-01-01")

# recoded data
print(tbl_df(filter(field.db1, field_subj_id == 12)))

# There was no subject 11, 28, 53, 73, or 76; Samples were generated in error.
# Delete samples for non-existant subjects
field.db1 <- filter(field.db1, field_subj_id != 11)%>%
  filter(field_subj_id != 28) %>%
  filter(field_subj_id != 53) %>%
  filter(field_subj_id != 73) %>% 
  filter(field_subj_id != 76)

# There was no second visit for subj 30; Samples generated in error.
# Delete samples for subject 30 on 2012-12-18. 

# Original data
print(tbl_df(filter(field.db1, field_subj_id == 30)))

# Make correction
field.db1 <- field.db1 %>%
  filter(!(field_subj_id == 30 & date_visit == as.Date("2012-12-18"))) 

# Corrected data subject 30
print(tbl_df(filter(field.db1, field_subj_id == 30)))

# There was no second visit for subj 120 on 2013-02-08 and sample 120_12 is not in REDCap sample database

# Original data
print(tbl_df(filter(field.db1, field_subj_id == 120)))

# Delete sample 120_12
field.db1 <- filter(field.db1, sample_id != "120_12")

# Corrected data
print(tbl_df(filter(field.db1, field_subj_id == 120)))

# Delete  samples (NP only) from erroneous second g2 visit for subject 69 (see above)
field.db1 <- filter(field.db1, sample_id != "69_6" & sample_id != "69_7")

## End of Data Editing Field Sample Database ##

# Number of columns in EDITED filed sample database (selected columns)
print(ncol(field.db1))

# Variable names in field.db1
print(head(tbl_df(field.db1)))

# How many subjectIDs are represented in the field.db1?
field.db1_subjectID_check <- field.db1 %>%
  group_by(field_subj_id) %>%
  count()
print(nrow(field.db1_subjectID_check))
# Data for 355 individuals - all samples collected - including clinic and gii

#### MERGE FIELD SAMPLE DATABASE WITH COMBINED CLINICAL DATABASE & G2 LOG ####

## Merge2 = merge of merge1 with field.db1 by field_subj_id and date_visit ##
merge2 <- merge(merge1, field.db1, by = c("field_subj_id", "date_visit"), all = T)

# Check number of subjectIDs now - should still be 355
merge2_subjectID_check <- merge2 %>%
  group_by(field_subj_id) %>%
  count()
print(nrow(merge2_subjectID_check))
# Yes, this seems to be correct. 

# Head of merge2
print(head(tbl_df(merge2)))

# Source of data in rows of merge2
print(ftable(addmargins(with(merge2, table(merge1.i, field.db1.i, exclude = c())))))

# Remove all samples that were assigned in the field db but not collected from the unenrolled subjects
merge2 <- merge2 %>% 
  filter(!(enrolled == F & sample_type %in% 
             c("GII condensate NO mask", "Throat Swab", "Impactor 5 um NO mask", "anterior nasal swab")))

# Source of data in rows of merge2 after removing extraneous samples
print(ftable(addmargins(with(merge2, table(merge1.i, field.db1.i, exclude = c())))))

# Merge2: rows where merge1 was not matched by rows from field.db1
print(filter(merge2, is.na(field.db1.i)))

# All rows in merge2 for subjects that have some merge1 rows not matching field.db1
x <- distinct(select(filter(merge2, is.na(field.db1.i)), field_subj_id))
print(inner_join(merge2, x, by = "field_subj_id"))

# All rows in field.db1 for subjects that had some merge1 rows not matching field.db1
print(inner_join(field.db1, x, by = "field_subj_id"))

# Subject 135 returned for a second screening visit but was never enrolled. 
# There is nothing in the REDCap clinical record to explain why no samples in the field DB were associated with the second visit.
# There are also no samples in the REDCap sample database for subject 135 except 135_1.

# As a result, we delete field_subj_id == 135 & date_visit == 2013-02-05 
merge2 <- merge2 %>%
  filter(!(field_subj_id == 135 & date_visit == "2013-02-05"))

# Checking merge2 rows where field.db1 not matched by rows from merge1
print(filter(merge2, is.na(merge1.i)))

# Checking: merge1 rows for subjects who had some field.db1 rows not matching merge1 rows
x <- distinct(select(filter(merge2, is.na(merge1.i)), field_subj_id))
nrow(x)

# Checking: Table of these rows
print(tbl_df(inner_join(merge1, x, by = "field_subj_id")))

# Source of data in rows of merge2 after removing 135_6.
print(ftable(addmargins(with(merge2, table(merge1.i, field.db1.i, exclude = c())))))

# Giving an indicator variable to this finalized merge2 df
merge2$merge2.i <- T

# Head of merge2
print(head(tbl_df(merge2)))

# Repeating from above now that data manipulation has been done: Check number of subjectIDs now - should still be 355
merge2_subjectID_check <- merge2 %>%
  group_by(field_subj_id) %>%
  count()
print(nrow(merge2_subjectID_check))
# Yes, this seems to be correct. 

#### READ in and work with the UMD SAMPLES DATABASE (REDCAP DATA) ####
sample_in_file <- 'UMD_Raw_Data/REDCAP/EMITUMDSamples2013_DATA.csv'
sample_in <- read.csv(sample_in_file, as.is = T)

sample_in$count_tech <- as.factor(sample_in$count_tech)

# Input UMD samples file (from REDCap)
sample_in_file

# Check number of subjectIDs now
sample_in_subjectID_check <- sample_in %>%
  group_by(field_subj_id) %>%
  filter(!is.na(field_subj_id)) %>%
  summarise(count = n())
print(nrow(sample_in_subjectID_check))
# Still data for all 355 screened participants.

# Number of rows
print(nrow(sample_in))
# Number of cols
print(ncol(sample_in))

sample_in <- sample_in %>%
  rename(date_visit = dt_visit)
sample_in$date_visit <- as.Date(sample_in$date_visit, format = "%m/%d/%Y")

# Head of sample_in dataframe after renaming date_visit and setting count_tech as factor
print(tbl_df(sample_in))

# Number of rows that have a missing sample_id
print(sum(sample_in$sample_id == ""))

# Number of rows that have a NA for sample_id
print(sum(is.na(sample_in$sample_id)))

## Data Editing ##
# Samples 69_6 and 69_7 collected for a g2_run = 2 that was not performed. See above data editing section. Deleted here.
# Also, samples 20-1 & 97-9 are typos and duplications. They are also deleted here.

sample_in <- sample_in %>%
  filter(!(sample_id %in% c("69_6", "69_7", "20-1", "97-9")))
sample_in$dt_stained <- with(sample_in, ifelse(sample_id == "20_1" & dt_count == "12/15/2012", "12/15/2012", dt_stained))

collection <- sample_in %>%
  filter(redcap_event_name == "collection_arm_1") %>%
  select(sample_id, field_subj_id, sample_type, date_visit)

assay1 <- sample_in %>% 
  filter(grepl('^assay1', redcap_event_name))

assay2 <- sample_in %>% 
  filter(grepl('^assay2', redcap_event_name))

# Number of collection records
print(nrow(collection))
# Number of assay 1 records
print(nrow(assay1))
# Number of assay 2 records
print(nrow(assay2))

# Event name and sample_type read in from REDCap sample database
print(ftable(addmargins(with(sample_in, table(redcap_event_name, sample_type, exclude = c())))))

# Head of sample_in samples with no sample_type
print(head(filter(sample_in, sample_type == "")))

passage <- select(filter(sample_in, !is.na(passage_1)), 
                  sample_id, 
                  passage_id, 
                  passage_id_problem, 
                  dt_pass, 
                  pass_tech, 
                  passage_1, 
                  passage_2, 
                  dt_pass_2, 
                  passage_complete)

## PROBLEM OBSERVATIONS THAT NEED EDITING ##

# Sample_id does not match Passage_id. Refer for lab review. 
# Meanwhile, will use sample_id as it is not duplicated.
print(tbl_df(filter(passage, sample_id != passage_id)))

passage$passpos <- (passage$passage_1 == 2 | passage$passage_2 == 2) # Passage is + if either passage is +
passage$validp <- !is.na(passage$passpos)
passage <- passage %>%
  select(sample_id, passpos, validp)

# Initial look at passage assays
# Number of passage assays
print(nrow(passage))

# Number of passage assays with missing date of passage
print(sum(passage$dt_pass == ""))

# Number of valid passage assays
print(sum(passage$validp))

# Number of invalid passage assays
print(sum(!passage$validp))

# Number of positive passage assays
print(sum(passage$passpos, na.rm = T))

# Number of negative passage assays
print(sum(!passage$passpos, na.rm = T))

focus1 <- filter(assay1, !dt_count == "")[ , c(1, 17:37)]

focus2 <- filter(assay2, !dt_count == "")[ , c(1, 17:37)]

# Focus Assays

# Samples with miss match sample_id and focus_id focus1
print(sum(!(focus1$sample_id == focus1$focus_id)))

# Samples with miss match sample_id and focus_id focus2
print(sum(!(focus2$sample_id == focus2$focus_id)))

# Number of focus1 counts
print(nrow(focus1))
# Number of focus2 counts
print(nrow(focus2))

# Number of focus1 rows with no dt_count
print(sum(focus1$dt_count == ""))
# Number of focus2 rows with no dt_count
print(sum(focus2$dt_count == ""))

#### Computation of focus assay results ####

focus1$df <- 10^(ifelse(is.na(focus1$dilution_factor), 0, focus1$dilution_factor))

focus2$df <- 10^(ifelse(is.na(focus2$dilution_factor), 0, focus2$dilution_factor))

area.24 <- pi*(15.4/2)^2
area.g <- 0.64

focus1$ct_24g <- rowSums(focus1[ , c(11:20)], na.rm = T) / (10*area.g)*area.24*focus1$df / 150*1000
focus1$ct_24w <- focus1$well*focus1$df / 150*1000
focus1$ct_96  <- rowSums(focus1[ , c(11:13)], na.rm = T)*focus1$df / 150*1000

focus1_24g <- focus1 %>%
  filter((focus1$plate_type == 1 | is.na(focus1$plate_type)) & focus1$count_meth == 1) %>% 
  select(-ct_96, -ct_24w) %>%
  rename(ct = ct_24g)

focus1_24w <- focus1 %>% 
  filter((focus1$plate_type == 1 | is.na(focus1$plate_type)) & focus1$count_meth == 2) %>% 
  select(-ct_96, -ct_24g) %>%
  rename(ct = ct_24w)

focus1_96  <- focus1 %>%
  filter(focus1$plate_type == 2) %>% 
  select(-ct_24w, -ct_24g) %>%
  rename(ct = ct_96)

focus1_c <- arrange(rbind(focus1_96, focus1_24w, focus1_24g))

focus2$ct_24g <- rowSums(focus2[ , c(11:20)], na.rm = T) / (10*area.g)*area.24*focus2$df / 150*1000

focus2$ct_24w <- focus2$well*focus2$df / 150*1000

focus2$ct_96  <- rowSums(focus2[ , c(11:13)], na.rm = T)*focus2$df / 150*1000

focus2_24g <- focus2 %>%
  filter((focus2$plate_type == 1 | is.na(focus2$plate_type)) & focus2$count_meth == 1) %>% 
  select(-ct_96, -ct_24w) %>%
  rename(ct = ct_24g)

focus2_24w <- focus2 %>%
  filter((focus2$plate_type == 1 | is.na(focus2$plate_type)) & focus2$count_meth == 2) %>% 
  select(-ct_96, -ct_24g) %>%
  rename(ct = ct_24w)

focus2_96  <- focus2 %>% 
  filter(focus2$plate_type == 2) %>% 
  select(-ct_24w, -ct_24g) %>%
  rename(ct = ct_96)

focus2_c <- arrange(rbind(focus2_96, focus2_24w, focus2_24g))

focus1_c <- select(focus1_c, sample_id, dt_count, count_tech, ct)

focus2_c <- select(focus2_c, sample_id, dt_count, count_tech, ct)

focus <- merge(focus1_c, focus2_c, by = "sample_id", all = T)
focus$ct <- rowMeans(cbind(focus$ct.x,focus$ct.y), na.rm = T)

missing_focus <- tbl_df(filter(focus, is.nan(ct)))
focus_allv <- focus
focus <- focus %>%
  select(sample_id, ct)

## FOCUS ASSAY RESULTS ##

# Samples listed as having a focus assay but without results
print(missing_focus)
summary(focus)

#### MERGE CULTURE RESULTS PIECE FROM FIELD SAMPLE DATABASE TO THE CUMULATIVE CLIN DB + G2 LOG + FIELD SAMPLE DB ####

# How many subjectIDs in the collection df?
collection_subjectID_check <- collection %>%
  group_by(field_subj_id) %>%
  filter(!is.na(field_subj_id)) %>%
  count()
print(nrow(collection_subjectID_check))
# All 355 are accounted for here

# How many subjectIDs in the passage df?
passage_subjectID_check <- passage %>%
  mutate(subject_id = gsub('_[0-9]*', '', sample_id)) %>%
  group_by(subject_id) %>%
  filter(!is.na(subject_id)) %>%
  count()
print(nrow(passage_subjectID_check))
# 158 are included in the passage data - not sure which 158 these are

# How many subjectIDs in the focus df?
focus_subjectID_check <- focus %>%
  mutate(subject_id = gsub('_[0-9]*', '', sample_id)) %>%
  group_by(subject_id) %>%
  filter(!is.na(subject_id)) %>%
  count()
print(nrow(focus_subjectID_check))
# 153 are included in the passage data - not sure which 153 these are

culture_results <- merge(collection, passage, by = "sample_id", all = T)
culture_results <- merge(culture_results, focus, by = "sample_id", all = T)

# N rows collection
print(nrow(collection))
# N rows passage
print(nrow(passage))
# N rows focus
print(nrow(focus))
# N rows culture_results
print(nrow(culture_results))

culture_results$np <- ifelse(culture_results$sample_type == 'Nasopharyngeal swab', T, F)
culture_results$impactor <- ifelse(culture_results$sample_type == 'Impactor 5 um NO mask', T, F)
culture_results$condensate <- ifelse(culture_results$sample_type == 'GII condensate NO mask', T, F)
culture_results$antnasal <- ifelse(culture_results$sample_type == 'anterior nasal swab', T, F)  
culture_results$throat <- ifelse(culture_results$sample_type == 'Throat Swab', T, F)  
culture_results$focus.i <- ifelse(is.na(culture_results$ct), F, T)
culture_results$passage.i <- ifelse(is.na(culture_results$validp), F, T)
culture_results$cr.i <- TRUE #Indicator for record present in culture_results

# Samples with (1) and without (0) passage assays and focus assays
print(ftable(addmargins(with(culture_results, table(sample_type, passage.i, focus.i, exclude = c()))))) 
## Note: All samples should have a sample type; anterior nasal swabs & impactors should not have culture assays ##

## PROBLEMATIC SAMPLES THAT NEED TO BE REVIEWED IN NOTEBOOKS AND REDCAP ##

# Based on review of culture results alone, may be resolved after merge with field field.db

# Samples with missing sample_type or sample_type=NA
print(tbl_df(filter(culture_results, is.na(sample_type)|sample_type == "")))

# Ant Nasal samples with either passage or focus assay results (even if neg/0)
print(tbl_df(culture_results %>% 
               filter(antnasal == T, focus.i == T | passage.i == T)))

# Impactor samples with either passage or focus assay results (even if neg/0)
print(tbl_df(culture_results %>% 
               filter(impactor == T, focus.i == T | passage.i == T)))

# NP swabs with a focus assay but no passage (even invalid) assay, or with a passage but no focus assay
print(tbl_df(culture_results %>% 
               filter(np == T, (focus.i == T & passage.i == F) | (focus.i == F & passage.i == T))))

# G-II condensate with a focus assay but no passage (even invalid) assay, or with a passage but no focus assay
print(tbl_df(culture_results %>% 
               filter(condensate == T, (focus.i == T & passage.i == F) | (focus.i == F & passage.i == T))))

# Culture_results
print(head(tbl_df((culture_results))))

## Merge3 = merge of merge2 with (culture_results from REDCap) by sample_id (only) ##
# x=merge2, y=culture_results

# Variables merge2
print(names(merge2))

# Variables culture_results
print(names(culture_results))

# How many subjectIDs in the collection df?
culture_results_subjectID_check <- culture_results %>%
  group_by(field_subj_id) %>%
  filter(!is.na(field_subj_id)) %>%
  count()
print(nrow(culture_results_subjectID_check))
# All 355 are accounted for here

merge3 <- merge(merge2, culture_results, c('sample_id'), all = TRUE)

# Number of rows in merge2
print(nrow(merge2))
# Number of rows in culture_results
print(nrow(culture_results))
# Number of rows in merge3
print(nrow(merge3))
#Number of cols in merge3
print(ncol(merge3))

print(head(tbl_df(merge3)))

# Table to check source of records after merge 3
print(ftable(addmargins(with(merge3, table(merge2.i, cr.i, exclude = c())))))

# Do date_visit match?
d.err <- filter(merge3, date_visit.x != date_visit.y | is.na(date_visit.x != date_visit.y))

# Number of samples where the date_visit.x (merge2) not equal date_visit.y (culture_results)
print(nrow(d.err))

# First 10 rows with non matching date_visit ordered by sample_id
print(top_n(select(d.err, 
                   sample_id, 
                   date_visit.x, 
                   date_visit.y, 
                   enrolled, 
                   visit_num, 
                   sample_type.x, 
                   sample_type.y), 
            10, sample_id))

# Do sample_type match?
t.err <- merge3 %>%
  filter(sample_type.x != sample_type.y | is.na(sample_type.x != sample_type.y))

# Number of rows where sample types don't match =
print(nrow(t.err))

# Columns show whether culture_result sample types were missing?, Rows likewise for merge2 sample type.
print(ftable(addmargins(with(merge3, 
                             table(miss.x <- is.na(sample_type.x), miss.y <- is.na(sample_type.y))))))

# Table of sample types by data source. (x=merge2, y=culture_results)
print(ftable(addmargins(with(merge3, 
                             table(merge2.i, cr.i, sample_type.x, sample_type.y, exclude = c())))))

# All non-matching sample types seem to be due to missing (NA) values.

# All sample ids begining with 237 in merge3
print(filter(merge3, grepl("^237", sample_id)))

# All sample ids begining with 237 in culture_results
print(filter(culture_results, grepl("^237", sample_id)))

# Sample 237_6 is an enrolled roommate, enrolled based on fever, therefore second NP swab should be in the lab.

# Samples in merge3 where culture_results data had no match in merge2 (field data)
x <- filter(merge3, cr.i == T, is.na(merge2.i))
print(select(x, sample_id))

# All sample ids begining with 301 in merge3
print(filter(merge3, grepl("^301", sample_id)))

# All records in merge3 with field_subj_id 301 from either source dataframe
print(filter(merge3, field_subj_id.x == 301 | field_subj_id.y == 301))

# Records for samples 301_x in culture_results: \n")
print(filter(culture_results, grepl("^301", sample_id)))

# All samples collected on same date as subject 301

x <- select(filter(merge3, field_subj_id.y == "301"), date_visit.x)

y <- right_join(merge3, x, by = "date_visit.x")

# Throat and Condensate samples from enrolled subject on the same day as 301
print(filter(y, sample_type.x %in% c("Throat Swab", "GII condensate NO mask")))

# Looks like 301_3 and 301_5 are erroneous duplicative entries for 302_3 and 302_5: Will delete extra 301 samples.
merge3 <- merge3 %>%
  filter(!(sample_id %in% c("301_3", "301_5")))

# Table to check source of records after merge 3 clean-up
print(ftable(addmargins(with(merge3, table(merge2.i, cr.i, exclude = c())))))

merge3 <- mutate(
  merge3, 
  subject_id = ifelse((!is.na(field_subj_id.x) | field_subj_id.x == ""), 
                      field_subj_id.x, 
                      field_subj_id.y),
  sample_type = ifelse((sample_type.x %in% 
                          c("Nasopharyngeal swab",
                            "Impactor 5 um NO mask", 
                            "GII condensate NO mask", 
                            "anterior nasal swab",
                            "Throat Swab")),
                       sample_type.x, 
                       sample_type.y),
  date_visit =  as.Date(ifelse(!is.na(date_visit.x), date_visit.x, date_visit.y), origin = "1970-01-01")
)

merge3 <- merge3 %>%
  select(-contains(".x"), -contains(".y"))

# Number of rows in merge3 after clean-up with no sample type
print(nrow(filter(merge3, is.na(sample_type) | sample_type == "")))

# Number of rows in merge3 after clean-up with no subject id
print(nrow(filter(merge3, is.na(subject_id) | subject_id == "")))

# Number of rows in merge3 after clean-up with no date visit
print(nrow(filter(merge3, is.na(date_visit))))

# Samples without sample type
print(filter(merge3, is.na(sample_type)))
# What to do with these?

## For the group not enrolled, figure out how to keep the one NP swabs that was sent to the lab and discard the other record.

# Examine the first visit to see which NP samples were cultured

np <- select(filter(merge3, np == T|sample_type == "Nasopharyngeal swab"), 
             subject_id, 
             sample_id, 
             enrolled, 
             visit_num, 
             validp, 
             ct) 

np1 <- np %>% 
  filter(visit_num == 1)

np1a <- select(mutate(np1, cultured = !is.na(validp) | !is.na(ct)), 
               sample_id,
               subject_id,
               cultured)

np1a <- select(separate(np1a, sample_id, c("id", "sample.id"), "_"), 
               subject_id, 
               sample.id, 
               cultured)

np1a.s <- spread(np1a, sample.id, cultured)

names(np1a.s)[2:length(names(np1a.s))] <- paste("sample", names(np1a.s)[2:length(names(np1a.s))], sep = "_")

print(summary(np1a.s))

# NA, means not cultured / assayed. 
# Conclude that all cultured visit 1 samples were called \"_1\" by the lab.

# Examine all visits for which NP samples were cultured

# Before spread can run on all np samples, must correct for reassigning 225 to 250 to avoid duplication of sample_#

np <- np %>%
  mutate(sample_id = 
           ifelse(subject_id == 250 & visit_num == 2, paste(sample_id, "a", sep = ""), sample_id)
         )

np.a <- select(mutate(np, cultured = !is.na(validp)|!is.na(ct)), sample_id, subject_id, cultured)
np.a <- select(separate(np.a, sample_id, c("id", "sample.id"), "_"), subject_id, sample.id, cultured)
np.a.s <- spread(np.a, sample.id, cultured)
names(np.a.s)[2:length(names(np.a.s))] <- paste("sample", names(np.a.s)[2:length(names(np.a.s))], sep = "_")
print(summary(np.a.s))

# Subjects with sample_6 cultured or subject_id = 250 (after reassigning subject 225 to subject 250)
print(filter(np.a.s, sample_6 == T|subject_id == 250))

# All samples for subject 247
print(filter(merge3, 
             subject_id == "247"))

# All samples for subje_id == 250
print(filter(merge3, 
             subject_id == "250"))

merge3 <- merge3 %>%
  rename( subject.id = subject_id, 
          sample.id = sample_id, 
          focus.ct = ct, 
          sample.type = sample_type, 
          date.visit = date_visit, 
          g2.run = g2_run, 
          visit.num = visit_num)

samples.cc <- merge3 %>%
  select(subject.id, 
         date.visit, 
         date_on_sx,
         sample.id, 
         sample.type, 
         g2.run, 
         visit.num, 
         passpos, 
         validp, 
         focus.ct, 
         enrolled, 
         rapid_flu___3, 
         rapid_flu_loc, 
         body_temp, 
         cough_number,
         sneeze_number)

# Check the number of subjectIDs in the samples.cc df
samples.cc_subjectID_check <- samples.cc %>%
  group_by(subject.id) %>%
  filter(!is.na(subject.id)) %>%
  count()
print(nrow(samples.cc_subjectID_check))
# There are 355 subjectIDs here = full set of screened participants

#### Write out EMIT_samples.cc.RDS file from merge of Clin DB + G2 Log + Field Sample DB ####

saveRDS(samples.cc, file = "Curated Data/Cleaned Data/EMIT_samples.cc.RDS")

#### **** Using Script: Jing's "interun calibration.R" **** ####

### 
## Original File Information

# Author: Jing Yan
# Date: August 9th, 2016
# Title: Interun calibration
#	Purpose: calibrate the Ct values based on the standard interrun calibrator	
# Input files:InputFiles_UMD
# InputFiles_UMD/PCR results/2016.08.05 1st visit NP swab FluA quant.csv
# InputFiles_UMD/PCR results/2016.08.05 1st visit NP swab FluB quant.csv
# InputFiles_UMD/PCR results/2016.08.08 GII and repeat NP swab FluA part Ia.csv
# InputFiles_UMD/PCR results/2016.08.08 GII and repeat NP swab FluA part Ib.csv
# InputFiles_UMD/PCR results/2016.08.05 GII and repeat NP swab FluA Part II.csv
# InputFiles_UMD/PCR results/2016.08.05 GII and repeat NP swabs FluA part III.csv
# InputFiles_UMD/PCR results/2016.08.08 GII and repeat NP swab FluB part Ia.csv
# InputFiles_UMD/PCR results/2016.08.08 GII and repeat NP swab FluB part Ib.csv
# InputFiles_UMD/PCR results/2016.06.24 GII and repeat NP swab FluB part II.csv
# Output files: R_output

#FluA standard curve:
#Y = -3.143*LOG(X) + 37.14
#FluB standard curve:
#Y = -3.167*LOG(X) + 33.75
#FluA standard curve 1st NP swab
#Y = -3.346*LOG(X) + 37.52
#FluB standard curve 1st NP swab
#Y = -3.297*LOG(X) + 34.45

###

# I'll note here that the first visit swabs and the rest of the pcr samples received different calibrations, ...
# ... in addition to the A's and B's getting different calibrations
# All of the cleaned files are written out into the Curated Data/Cleaned Data directory

low = 53.2
high = 53200
ctLAII = -3.143*log10(low) + 37.14
ctHAII = -3.143*log10(high) + 37.14
ctLA1np = -3.346*log10(low) + 37.52
ctHA1np = -3.346*log10(high) + 37.52
ctLBII = -3.167*log10(low) + 33.75
ctHBII = -3.167*log10(high) + 33.75
ctLB1np = -3.297*log10(low) + 34.45
ctHB1np = -3.297*log10(high) + 34.45

#### READ in *"2016.08.05 1st visit NP swab FluA quant.csv"* ####
fluAnp <- read.csv('UMD_Raw_Data/PCR Data/PCR results/2016.08.05 1st visit NP swab FluA quant.csv', as.is = T)
fluAnp <- fluAnp %>% 
  filter(Ct..dRn. != 'Reference')
fluAnp$Ct..dRn. <- as.numeric(fluAnp$Ct..dRn.)

fluAnp1 <- fluAnp %>% 
  select(Experiment, Well.Name, Ct..dRn.)

fluAnplow <- fluAnp1 %>% 
  filter(grepl('_Low', Well.Name))
fluAnplow$dif <- fluAnplow$Ct..dRn.-ctLA1np

fluAnphigh <- fluAnp1 %>%
  filter(grepl('_High', Well.Name))
fluAnphigh$dif <- fluAnphigh$Ct..dRn.-ctHA1np

fluAnp2 <- rbind(fluAnplow, fluAnphigh) %>% 
  ungroup %>%
  arrange(Experiment) %>%
  mutate(date = gsub('^[0-9]*.', '', Experiment))

fluAnp3 <- fluAnp2 %>%
  group_by(date) %>%
  mutate(avgdiff = mean(dif))
fluAnp3$cfactor <- 10^((fluAnp3$avgdiff) / 3.346)
fluAnp3 <- fluAnp3 %>% 
  distinct(date, cfactor, .keep_all = TRUE)

fluAnp4 <- fluAnp3 %>% 
  select(date, cfactor)

saveRDS(fluAnp4, "Curated Data/Cleaned Data/fluA_1np_calibration.RDS")

#### READ in *"2016.08.05 1st visit NP swab FluB quant.csv"* ####
fluBnp <- read.csv('UMD_Raw_Data/PCR Data/PCR results/2016.08.05 1st visit NP swab FluB quant.csv', as.is = T)
fluBnp <- fluBnp %>%
  filter(Ct..dRn. != 'Reference')
fluBnp$Ct..dRn. <- as.numeric(fluBnp$Ct..dRn.)

fluBnp1 <- fluBnp %>% 
  select(Experiment, Well.Name, Ct..dRn.)

fluBnplow <- fluBnp1 %>% 
  filter(grepl('_Low', Well.Name))
fluBnplow$dif <- fluBnplow$Ct..dRn.-ctLB1np
fluBnplow <- fluBnplow %>%
  filter(!is.na(Ct..dRn.))

fluBnphigh <- fluBnp1 %>% 
  filter(grepl('_High', Well.Name))
fluBnphigh$dif <- fluBnphigh$Ct..dRn.-ctHB1np
fluBnphigh <- fluBnphigh %>% 
  filter(!is.na(Ct..dRn.))

fluBnp2 <- rbind(fluBnplow, fluBnphigh) %>% 
  ungroup %>%
  arrange(Experiment) %>%
  mutate(date = gsub('^[0-9]*.', '', Experiment))

fluBnp3 <- fluBnp2 %>% 
  group_by(date) %>%
  mutate(avgdiff = mean(dif))
fluBnp3$cfactor <- 10^((fluBnp3$avgdiff) / 3.297)
fluBnp3 <- fluBnp3 %>% 
  distinct(date, cfactor, .keep_all = TRUE)

fluBnp4 <- fluBnp3 %>%
  select(date, cfactor)

saveRDS(fluBnp4, "Curated Data/Cleaned Data/fluB_1np_calibration.RDS")

#### READ in and work with the 4 PCR raw datafiles for flu A (not including the '1st visit' file) ####
# READ in "2016.08.08 GII and repeat NP swab FluA part Ia.csv" 
# READ in "2016.08.08 GII and repeat NP swab FluA part Ib.csv" 
# READ in "2016.08.05 GII and repeat NP swab FluA Part II.csv" 
# READ in "2016.08.05 GII and repeat NP swabs FluA part III.csv" 

fluAI1 <- read.csv('UMD_Raw_Data/PCR Data/PCR results/2016.08.08 GII and repeat NP swab FluA part Ia.csv', as.is = T)
fluAI2 <- read.csv('UMD_Raw_Data/PCR Data/PCR results/2016.08.08 GII and repeat NP swab FluA part Ib.csv', as.is = T)
fluAII <- read.csv('UMD_Raw_Data/PCR Data/PCR results/2016.08.05 GII and repeat NP swab FluA Part II.csv', as.is = T)
fluAIII <- read.csv('UMD_Raw_Data/PCR Data/PCR results/2016.08.05 GII and repeat NP swabs FluA part III.csv', as.is = T)

fluAI <- rbind(fluAI1, fluAI2) %>% 
  ungroup()
fluAI <- fluAI %>% 
  filter(Ct..dRn. != 'Reference')
fluAI$Ct..dRn. <- as.numeric(fluAI$Ct..dRn.)
fluAI <- fluAI %>% 
  select(Experiment, Well.Name, Ct..dRn.)

fluAII <- fluAII %>% 
  filter(Ct..dRn. != 'Reference')
fluAII$Ct..dRn. <- as.numeric(fluAII$Ct..dRn.)

fluAII1 <- fluAII %>% 
  select(Experiment, Well.Name, Ct..dRn.)

fluAIII <- fluAIII %>% 
  filter(Ct..dRn. != 'Reference')
fluAIII$Ct..dRn. <- as.numeric(fluAIII$Ct..dRn.)

fluAIII1 <- fluAIII %>% 
  select(Experiment, Well.Name, Ct..dRn.)

fluA1 <- rbind(fluAI, fluAII1, fluAIII1) %>%
  ungroup()

fluA1low <- fluA1 %>%
  filter(grepl('_Low', Well.Name))
fluA1low$dif <- fluA1low$Ct..dRn.-ctLAII

fluA1high <- fluA1 %>%
  filter(grepl('_High', Well.Name))
fluA1high$dif <- fluA1high$Ct..dRn.-ctHAII

fluA2 <- rbind(fluA1low, fluA1high) %>% 
  ungroup()
fluA2 <- fluA2[order(fluA2$Experiment), ]
fluA2 <- fluA2 %>%
  mutate(date = gsub('^[0-9]*.', '', Experiment))

fluA3 <- fluA2 %>% 
  group_by(date) %>% 
  mutate(avgdiff = mean(dif))
fluA3$cfactor <- 10^((fluA3$avgdiff) / 3.143)
fluA3 <- fluA3 %>%
  distinct(date, cfactor, .keep_all = TRUE)

fluA4 <- fluA3 %>%
  select(date, cfactor)

saveRDS(fluA4, "Curated Data/Cleaned Data/fluA_calibration.RDS")

#### READ in and work with the 3 PCR raw datafiles for flu B (not including the '1st visit' file) ####
# Read in: "2016.08.08 GII and repeat NP swab FluB part Ia.csv"
# Read in: "2016.08.08 GII and repeat NP swab FluB part Ib.csv"
# Read in: "2016.06.24 GII and repeat NP swab FluB part II.csv"

fluBI1 <- read.csv('UMD_Raw_Data/PCR Data/PCR results/2016.08.08 GII and repeat NP swab FluB part Ia.csv', as.is = T)
fluBI2 <- read.csv('UMD_Raw_Data/PCR Data/PCR results/2016.08.08 GII and repeat NP swab FluB part Ib.csv', as.is = T)
fluBII <- read.csv('UMD_Raw_Data/PCR Data/PCR results/2016.06.24 GII and repeat NP swab FluB part II.csv', as.is = T)

fluBI <- rbind(fluBI1, fluBI2) %>% 
  ungroup
fluBI <- fluBI %>% 
  filter(Ct..dRn. != 'Reference')
fluBI$Ct..dRn. <- as.numeric(fluBI$Ct..dRn.)

fluBI1 <- fluBI %>%
  select(Experiment, Well.Name, Ct..dRn.)

fluBII <- fluBII %>%
  filter(Ct..dRn. != 'Reference')
fluBII$Ct..dRn. <- as.numeric(fluBII$Ct..dRn.)

fluBII1 <- fluBII %>%
  select(Experiment, Well.Name, Ct..dRn.)

fluB1 <- rbind(fluBI1, fluBII1) %>%
  ungroup

fluB1low <- fluB1 %>% 
  filter(grepl('_Low', Well.Name))
fluB1low$dif <- fluB1low$Ct..dRn.-ctLBII

fluB1high <- fluB1 %>%
  filter(grepl('_High', Well.Name))
fluB1high$dif <- fluB1high$Ct..dRn.-ctHBII

fluB2 <- rbind(fluB1low, fluB1high) %>% 
  ungroup
fluB2 <- fluB2[order(fluB2$Experiment),]
fluB2 <-fluB2 %>% 
  filter(!is.na(Ct..dRn.))
fluB2 <- fluB2 %>%
  mutate(date = gsub('^[0-9]*.','',Experiment))

fluB3 <- fluB2 %>%
  group_by(date) %>%
  mutate(avgdiff = mean(dif))
fluB3$cfactor <- 10^((fluB3$avgdiff) / 3.167)
fluB3 <- fluB3 %>%
  distinct(date, cfactor, .keep_all = TRUE)

fluB4 <- fluB3 %>%
  select(date, cfactor)

saveRDS(fluB4, "Curated Data/Cleaned Data/fluB_calibration.RDS")

#### **** Using Script: Jing and Dr. Milton's "np fine coarse pcr quantity.R"**** ####

###
# Original file information:

# Author: Jing Yan& Don Milton
# Date: September 17, 2015
# Revision Date: Jan 24, 2016
# Title: np fine coarse pcr quantity.R
# Purpose: To sort the data files from the lab (Michael Grantham) PCR data for 2rd or 3rd np, fine and coarse 

# Input files:InputFiles_UMD
# InputFiles_UMD/PCR_1.8.2016/2016.1.4 GII and repeat NP samples FluA Part I.csv
# InputFiles_UMD/PCR_1.8.2016/2016.1.4 GII and repeat NP samples FluA Part II.csv
# InputFiles_UMD/PCR_1.8.2016/2016.1.8 GII samples and repeat NP swabs FluB Part I.csv
# InputFiles_UMD/PCR_1.8.2016/2016.1.8 GII samples and repeat NP swabs FluB Part II.csv

# Output files: R_output

# Question: How to treat the samples with multiple PCR results? 

###

#### READ in "fluA_calibration.RDS" ####
fluAcali <- readRDS("Curated Data/Cleaned Data/fluA_calibration.RDS")

#### READ in "2016.08.08 GII and repeat NP swab FluA part Ia.csv" ####
pcr_A1a <- read.csv('UMD_Raw_Data/PCR data/PCR results/2016.08.08 GII and repeat NP swab FluA part Ia.csv', as.is = T)

#### READ in "2016.08.08 GII and repeat NP swab FluA part Ib.csv" ####
pcr_A1b <- read.csv('UMD_Raw_Data/PCR data/PCR results/2016.08.08 GII and repeat NP swab FluA part Ib.csv', as.is = T)

# Bind the two flu A files together
pcr_A1 <- rbind(pcr_A1a, pcr_A1b) %>%
  ungroup()

names(pcr_A1)

# Clean the flu A data
pcr_A1 <- pcr_A1 %>%
  filter(Ct..dRn. != 'Reference') %>%
  filter(!grepl('_Low', Well.Name)) %>%
  filter(!grepl('NTC', Well.Name)) %>%
  filter(!grepl('Standard', Well.Type)) %>%
  filter(!grepl('High', Well.Name)) %>%
  filter(grepl('_', Well.Name)) %>%
  select(-Well, -Well.Type, -Threshold..dRn., -Replicate) %>%
  mutate(Well.Name = gsub('_A','', Well.Name)) %>%
  mutate(Well.Name = gsub('A_','', Well.Name)) %>%
  mutate(subject.id = gsub('_[0-9]*', '', Well.Name)) %>%
  rename(copies.in = Quantity..copies.)
pcr_A1[pcr_A1 == "No Ct"] <- ''
pcr_A1$copies.in <- as.numeric(pcr_A1$copies.in)

# Number of rows in part1 (gii and 2rd or 3rd NP samples) influenza A PCR data
nrow(pcr_A1)
# Number of columns in part1 (gii and 2rd or 3rd NP samples) influenza A PCR data
ncol(pcr_A1)

pcr_A1 <- pcr_A1[order(pcr_A1$Well.Name), ]

#### READ in "2016.08.05 GII and repeat NP swab FluA Part II.csv" ####

pcr_A2 <- read.csv('UMD_Raw_Data/PCR data/PCR results/2016.08.05 GII and repeat NP swab FluA Part II.csv',as.is=T)
names(pcr_A2)

pcr_A2 <- pcr_A2 %>%
  filter(Ct..dRn. != 'Reference') %>%
  filter(!grepl('_Low', Well.Name)) %>%
  filter(!grepl('NTC', Well.Name)) %>%
  filter(grepl('_', Well.Name)) %>%
  filter(!grepl('High', Well.Name)) %>%
  filter(!grepl('Standard', Well.Type)) %>%
  mutate(Well.Name = gsub('_A', '', Well.Name)) %>%
  mutate(Well.Name = gsub('A_', '', Well.Name))%>%
  rename(copies.in = Quantity..copies.) %>%
  mutate(subject.id = gsub('_[0-9]*', '', Well.Name)) %>% 
  select(-Well, -Well.Type, -Threshold..dRn., -Replicate, -X, -FAM..Y....3.143.LOG.X....37.14..Eff....108.0.)
pcr_A2[pcr_A2 == "No Ct"] <- ''
pcr_A2$copies.in <- as.numeric(pcr_A2$copies.in)

# Number of rows in part2 (gii and 2rd or 3rd NP samples) influenza A PCR data
nrow(pcr_A2)
# Number of columns in part2 (gii and 2rd or 3rd NP samples) influenza A PCR data
ncol(pcr_A2)

pcr_A2 <- pcr_A2[order(pcr_A2$Well.Name), ]

#### READ in "2016.08.05 GII and repeat NP swabs FluA part III.csv" ####

pcr_A3 <- read.csv('UMD_Raw_Data/PCR data/PCR results/2016.08.05 GII and repeat NP swabs FluA part III.csv', as.is = T)
names(pcr_A3)
pcr_A3 <- pcr_A3 %>% 
  filter(Ct..dRn. != 'Reference') %>% 
  filter(!grepl('_Low', Well.Name)) %>%
  filter(!grepl('NTC', Well.Name)) %>%
  filter(grepl('_', Well.Name)) %>%
  filter(!grepl('Standard', Well.Type)) %>%
  filter(!grepl('B', Well.Name)) %>%
  filter(!grepl('High', Well.Name)) %>%
  mutate(Well.Name = gsub('_A', '', Well.Name)) %>%
  mutate(Well.Name = gsub('A_', '', Well.Name)) %>%
  rename(copies.in = Quantity..copies.) %>%
  select(-Well, -Well.Type, -Threshold..dRn.) %>%
  mutate(subject.id = gsub('_[0-9]*', '', Well.Name))
pcr_A3[pcr_A3 == "No Ct"] <- ''
pcr_A3$copies.in <- as.numeric(pcr_A3$copies.in)

# Number of rows in part2 (gii and 2rd or 3rd NP samples) influenza A PCR data
nrow(pcr_A3)
# Number of columns in part2 (gii and 2rd or 3rd NP samples) influenza A PCR data
ncol(pcr_A3)

pcr_A3 <- pcr_A3[order(pcr_A3$Well.Name), ]

#### Bind the flu A data ####

pcr_A <- rbind(pcr_A1, pcr_A2, pcr_A3) %>% 
  ungroup()
pcr_A <- pcr_A %>% 
  mutate(type = 'A') %>% 
  mutate(copy.num = copies.in*250)

# Number of rows of total (gii and 2rd or 3rd NP samples) influenza A PCR data
nrow(pcr_A)
# Number of columns of total (gii and 2rd or 3rd NP samples) influenza A PCR data
ncol(pcr_A)

pcr_A$subject.id <- as.numeric(pcr_A$subject.id)

pcr_A <- pcr_A %>%
  arrange(subject.id) %>%
  mutate(date = gsub('^[0-9]*.', '', Experiment))

pcr_Afinal <- pcr_A %>%
  left_join(fluAcali, by = 'date')
pcr_Afinal$virus.copies <- pcr_Afinal$copy.num*pcr_Afinal$cfactor

#### READ in "fluB_calibration.RDS" ####
fluBcali <- readRDS("Curated Data/Cleaned Data/fluB_calibration.RDS")

#### READ in "2016.08.08 GII and repeat NP swab FluB part Ia.csv" ####
pcr_B1a <- read.csv('UMD_Raw_Data/PCR data/PCR results/2016.08.08 GII and repeat NP swab FluB part Ia.csv', as.is = T)

#### READ in "2016.08.08 GII and repeat NP swab FluB part Ib.csv" ####
pcr_B1b <- read.csv('UMD_Raw_Data/PCR data/PCR results/2016.08.08 GII and repeat NP swab FluB part Ib.csv', as.is = T)

#### Work with the fluB pcr data ####

pcr_B1 <- rbind(pcr_B1a, pcr_B1b) %>% 
  ungroup
names(pcr_B1)

## DATA EDITING ##
# It appears that the 20.2015.08.05.013.248 2012-13 samples FluA_B and FluA PCR experiment is already included elsewhere as 17.2015.8.5.013.248 2012-13 samples FluA_B and FluA PCR. We will keep the 17.2015.8.5.013.248 2012-13 samples FluA_B and FluA PCR assay and eliminate the other.
# Also, it appears that the 15. 2015.07.29.013.243 2012-13 samples FluA_B and FluA PCR and 17. 2015.7.29.013.243 2012-13 samples FluA_B and FluA PCR assays are the same with just a minor change in the name. We will use the first of these, and eliminate the other
# Also, it appears that the 16. 2015.07.30.013.245 2012-13 samples FluA_B and FluA PCR and the 18. 2015.7.30.013.245 2012-13 samples FluA_B and FluA PCR dfs are carbon copies with minor naming changes. Let's elminat the second assay here.
pcr_B1 <- pcr_B1 %>%
  filter(Experiment != "20. 2015.08.05.013.248 2012-13 samples FluA_B and FluA PCR") %>%
  filter(Experiment != "17. 2015.7.29.013.243 2012-13 samples FluA_B and FluA PCR") %>%
  filter(Experiment != "18. 2015.7.30.013.245 2012-13 samples FluA_B and FluA PCR")

pcr_B1 <- pcr_B1 %>% 
  filter(Ct..dRn. != 'Reference') %>% 
  filter(!grepl('_Low', Well.Name)) %>% 
  filter(!grepl('NTC', Well.Name)) %>% 
  filter(grepl('_', Well.Name)) %>% 
  filter(!grepl('Standard', Well.Type)) %>% 
  filter(!grepl('High', Well.Name)) %>% 
  mutate(Well.Name = gsub('_B', '', Well.Name)) %>% 
  mutate(Well.Name = gsub('B_', '', Well.Name)) %>% 
  select(-Well.Type,-Replicate,-Threshold..dRn.,-Well) %>% 
  mutate(subject.id = gsub('_[0-9]*','',Well.Name)) %>%
  rename(copies.in = Quantity..copies.)
pcr_B1[pcr_B1 == "No Ct"] <- ''
pcr_B1$copies.in <- as.numeric(pcr_B1$copies.in)

# Number of rows in part1 (gii and 2rd or 3rd NP samples) influenza B PCR data
nrow(pcr_B1)
# Number of columns in part1 (gii and 2rd or 3rd NP samples) influenza B PCR data
ncol(pcr_B1)

#### READ in "2016.06.24 GII and repeat NP swab FluB part II.csv" ####
pcr_B2 <- read.csv('UMD_Raw_Data/PCR data/PCR results/2016.06.24 GII and repeat NP swab FluB part II.csv',as.is=T)
names(pcr_B2)

pcr_B2 <- pcr_B2 %>%
  filter(Ct..dRn. != 'Reference') %>%
  filter(!grepl('_Low', Well.Name)) %>% 
  filter(!grepl('NTC', Well.Name)) %>% 
  filter(grepl('_', Well.Name)) %>% 
  filter(!grepl('Standard', Well.Type)) %>%
  filter(!grepl('High', Well.Name)) %>% 
  mutate(Well.Name = gsub('_B', '', Well.Name)) %>%
  mutate(Well.Name = gsub('B_', '', Well.Name)) %>%
  mutate(subject.id = gsub('_[0-9]*', '', Well.Name)) %>% 
  select(-Well.Type, -Threshold..dRn., -Well) %>%
  rename(copies.in = Quantity..copies.)
pcr_B2[pcr_B2 == "No Ct"] <- ''
pcr_B2$copies.in <- as.numeric(pcr_B2$copies.in)

# Number of rows in part2 (gii and 2rd or 3rd NP samples) influenza B PCR data
nrow(pcr_B2)
# Number of columns in part2 (gii and 2rd or 3rd NP samples) influenza B PCR data
ncol(pcr_B2)

#### Merge the flu B pcr data ####

pcr_B <- arrange(rbind(pcr_B1, pcr_B2))
pcr_B <- pcr_B %>%
  mutate(type = 'B')
pcr_B <- pcr_B[order(pcr_B$Well.Name), ]
pcr_B <- pcr_B %>%
  mutate(copy.num = copies.in*272)

# Number of rows of total (gii and 2rd or 3rd NP samples) influenza B PCR data
nrow(pcr_B)
# Number of columns of total (gii and 2rd or 3rd NP sampless) influenza B PCR data
ncol(pcr_B)

pcr_B <- pcr_B[order(pcr_B$subject.id), ]
pcr_B <- pcr_B %>% 
  mutate(date = gsub('^[0-9]*.', '', Experiment))

pcr_Bfinal <- pcr_B %>%
  left_join(fluBcali, by = 'date')
pcr_Bfinal$virus.copies <- pcr_Bfinal$copy.num*pcr_Bfinal$cfactor

#### Merge the fluA and fluB pcr data ####

total.pcr <- rbind(pcr_Afinal, pcr_Bfinal)
total.pcr$Well.Name[total.pcr$Well.Name == '42_16'] <- '42_18'
total.pcr$Well.Name[total.pcr$Well.Name == '2110_2'] <- '210_2'
total.pcr$subject.id[total.pcr$subject.id == '2110'] <- '210'
total.pcr$Well.Name[total.pcr$Well.Name == '297_10'] <- '297_9'
total.pcr$Well.Name[total.pcr$Well.Name == '194_10'] <- '194_3'
total.pcr$Well.Name[total.pcr$Well.Name == '27_12'] <- '27_11'
total.pcr$Well.Name[total.pcr$Well.Name == '117_4'] <- '117_3'
total.pcr$Well.Name[total.pcr$Well.Name == '161_8'] <- '161_11'
total.pcr$Well.Name[total.pcr$Well.Name == '292_11' & total.pcr$Ct..dRn. == 30.82] <- '292_9'
total.pcr$Well.Name[total.pcr$Well.Name == '292_11' & total.pcr$Ct..dRn. == 30.71] <- '292_9'

# There are a bunch of UK aerosol samples in this pcr file. They have "sd" in their Well.Name.
# I can remove these. 
total.pcr <- total.pcr %>% 
  filter(!grepl('sd', Well.Name))

# Convert the subject.id to numeric class and order from least to highest subject ID
total.pcr$subject.id <- as.numeric(total.pcr$subject.id)
total.pcr <- total.pcr %>%
  arrange(subject.id)

# Number of rows of total (gii and 2rd or 3rd NP samples)  PCR data
nrow(total.pcr)
# Number of columns of total (gii and 2rd or 3rd NP samples) PCR data
ncol(total.pcr)

#### Before merging PCR data with subtype, need the subtypes ####
## In order to do this, we will insert the subtype script

#### **** Using Script: Jing and Dr. Milton's "Subtype analysis.R" **** ####

###
## Original file information:

# "Subtype analysis.R
# by Jing Yan & Don Milton
# Jan 24, 2016 -  
# Purpose: Sorting and analyze the sample subtypes from Lab data. 
# ... The types were based on the subtypes table on box in ...
# ... folder: EMIT_Date_Analysis. 
#________________________________________________________
#  Procedures: 
#  1. We input the subtype part I and subtype part II from the folder InputFiles_UMD/PCR_9.16.2015 
#     Files: '2015.09.15 Subtyping part I.csv' and '2015.09.15 Subtyping part II.csv'
#  2. We combine the part1 and part2, the rows equal to the sum of part1 and part2
#  3. We remove the rows within the data (e.g. reference, calibration, etc.) not needed for analysis
#  4. Assign a new column for subject ID, and another column for sample type 
#  5. Assign any subtype with "No Ct" as False, and the other ones with value as True (after checking for very high Ct values)
#  6. Determine the sample type based on the combination of T and F
#  7. Assign each sample type a different integer (1-10),check if all the subject has a number assigned based on subtype
#  8. We remove the subject with duplicated data (with same subject ID and sample type)
#  9. For the subject with multiple experiments but different sample type, pick them out and do case by case analysis
# 10. Combine the selected experiments with the other subjects with unique experiment( final_subtype)
# 11. Output (save) dataframe as EMIT_subtypes.RDS with limited number of variables. 
#________________________________________________________

###

#### READ in and work with "2016.06.17 1st visit NP swab subtyping.csv" ####

part1 <- read.csv('UMD_Raw_Data/PCR Data/PCR Results/2016.06.17 1st visit NP swab subtyping.csv', as.is = T)
part1 <- part1 %>% 
  select(-Well, -Well.Type, -Threshold..dR.)

# Number of rows in part11 file
print(nrow(part1))
# Number of columns in part11 file
print(ncol(part1))

# Remove the rows within the data (e.g. reference, calibration, etc.) not needed for analysis
m1 <- part1 %>% 
  filter(Ct..dR. != 'Reference') %>% 
  filter(!grepl('[A-Z]_[A-Z]', Well.Name)) %>% 
  filter(!grepl('PosControl', Well.Name)) %>% 
  mutate(subject_id = gsub('_1_[0-Z]*[0-Z]*[0-Z]*', '', Well.Name)) %>% # gen new col for subject_ID which can be obtained from sample_ID 
  mutate(subject_id = gsub('_[0-Z]*', '', subject_id))
m1$Experiment <- as.factor(m1$Experiment)

# Number of subj-exp in merged data
print(nrow(m1 %>% 
             distinct(subject_id)))

# Identify each subtype e.g. A, B, H3,give a new column called types
m7 <- m1 %>% 
  mutate(types = gsub('[0-9]*_1_', '', Well.Name)) %>% 
  mutate(types = gsub('1_', '', types))

m7$subject_id <- as.numeric(m7$subject_id)

m7 <- m7 %>%
  arrange(subject_id) # order the data based on subject_id

# Number of rows in the sorted subtype file
print(nrow(m7))
# Number of columns in the sorted subtype file
print(ncol(m7))

# Examine the Ct Values by type
m7$Ct <- as.double(m7$Ct..dR.)
print(filter(m7, !is.na(Ct)) %>% 
        group_by(types) %>% 
        summarise(avg = mean(Ct), mx = max(Ct)))

# Reactions with very high Ct values (>= 40) -- double check these, if there are any:
suspects <- filter(m7, Ct >= 40)
print(suspects)

m7$pos <- with(m7, ifelse(Ct..dR. == "No Ct", F, T)) # Considered Positive if any Ct value

m7 <- m7 %>%
  arrange(subject_id, types, Ct)

m7.a <- m7 %>% 
  group_by(subject_id, types) %>% 
  summarise(n = n())

m7.b <- m7.a %>% 
  filter(n >= 2) # find subj with more than one result

m7.c <- m7.b %>%
  inner_join(m7, by = c('subject_id', 'types')) #get type data for the subj with more than one type result

## Getting a df with true or false listed for each subtype for each subject_id, where the vars are the subtypes and the observations are the subject_ids

m7.c <- m7.c %>%
  arrange(subject_id, types, Ct) %>% 
  select(-n) %>% 
  distinct(subject_id, types, .keep_all = TRUE)

m7.c2 <- m7 %>%
  anti_join(m7.c, by = c('subject_id', 'types'))

m7 <- m7.c %>%
  full_join(m7.c2) %>%
  arrange(subject_id, types, Ct)

m7$Experiment <- as.factor(m7$Experiment)

m7 <- m7 %>% 
  distinct(subject_id, types, .keep_all = T)

m7_1 <- m7 %>% 
  select(subject_id, types, pos) 

m7.s <- spread(m7_1, key = types, value = pos) # cols = types, expt = rows

# Number of subject - experiments
print(nrow(m7.s))
# This m7.s gives the final set of subject_ids with their subtype classification

## Adding some additional variables to the m7.s df to further classify the observations

m7.s$type.sub.H1 <- with(m7.s, ifelse(A & H1, T, F)) # if A and H1 then subtype = H1

m7.s$type.sub.H3 <- with(m7.s, ifelse(A & H3, T, F)) # if A and H3 then subtype = H3

m7.s$type.sub.PH1 <- with(m7.s, ifelse(A & PA & PH1, T, F))

m7.s$type.B <- with(m7.s, ifelse(B, T, F))

m7.s$type.H3N2.and.B <- with(m7.s, ifelse(B & type.sub.H3, T, F))

m7.s$type.H3N2.and.PH1 <- with(m7.s, ifelse(type.sub.H3 & type.sub.PH1, T, F))

m7.s$type.B.and.PH1 <- with(m7.s, ifelse(type.B & type.sub.PH1, T, F))

m7.s$type.sub.indet <- with(m7.s, ifelse(((A | H3 | H1 | PH1 | PA) # any one of the A reactions
                                          & !(type.sub.H3|type.sub.H1|type.sub.PH1)), T, F)) # & did not meet a subt def

m7.s$type.neg <- with(m7.s, ifelse(RP & # RP positive w/o another pos is a true negative reaction.
                                     !(A|B|H3|PH1|PA), T, F))

m7.s$type.badass <- with(m7.s, ifelse(!RP & # RP negative and everything else neg is a bad assay.
                                        !(A | H3 | H1 | PH1 | PA  | B ), T, F))

## Assign a number code to each possible combination of results (among those observed)

m7.s$num <- NA
m7.s$num[m7.s$type.sub.H1 == T] <- 1
m7.s$num[m7.s$type.sub.H3 == T] <- 2
m7.s$num[m7.s$type.sub.PH1 == T] <- 3
m7.s$num[m7.s$type.B == T] <- 4
m7.s$num[m7.s$type.neg == T] <- 5
m7.s$num[m7.s$type.H3N2.and.B == T] <- 6
m7.s$num[m7.s$type.H3N2.and.PH1 == T] <- 7
m7.s$num[m7.s$type.B.and.PH1 == T] <- 8
m7.s$num[m7.s$type.sub.indet == T] <- 9
m7.s$num[m7.s$type.badass == T] <- 10

# Number of rows without a number assigned
print(nrow(filter(m7.s, is.na(num))))
# Number of rows without a number assigned - assigned to an object
check <- m7.s %>% 
  filter(is.na(num))

# Create a subset of the data (m7.s1) with subj who have only one assay or have more than one assay but the same result each time having only one row. And, with more than one row for those subj with different results for repeat assays
# Keep only rows that are different in both result 'num' and subject name
m7.s1 <- m7.s[order(m7.s$subject_id), ]

# Count the number of rows n for each subject_id and find out the subjects with multiple experiments but different sample type
m7.s2 <- m7.s1 %>% 
  group_by(subject_id) %>% 
  summarise(n = n())
# Number of subjects with 1, 2, or more obs that are different additional obs for single subjects that were not different have been deleted
print(with(m7.s2, addmargins(table(n, exclude = c()))))

## Move to classify a final subtype for each subject_id

finalsubtype <- m7.s1 %>%
  rename(subject.id = subject_id, type.inf = num) %>% 
  select(subject.id, type.inf, A, B, H1, H3, PH1, PA, RP)
finalsubtype$subject.id <- as.integer(finalsubtype$subject.id)

# Assign labels to type.inf
finalsubtype$type.inf <- 
  factor(finalsubtype$type.inf, levels = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10),
         labels = c('seasonal H1','H3N2','Pandemic H1','B','Negative','H3N2 and B','H3N2 and PH1','B and PH1','Indeterminate','bad assay'))
finalsubtype$type.inf <- as.character(finalsubtype$type.inf)

indeterminate <- finalsubtype %>% 
  filter(type.inf == 'Indeterminate')

finalsubtype$type.inf[finalsubtype$subject.id == 95] <- 'B and unsubtypable A'
finalsubtype$type.inf[finalsubtype$subject.id == 176] <- 'H3N2' 
finalsubtype$type.inf[finalsubtype$subject.id == 335] <- 'B' 
finalsubtype$type.inf[finalsubtype$subject.id == 64] <- 'Unsubtypable A' 

# There are 183 entries in the finalsubtype df
# Let's check how many negative, indeterminate, and others there are - also, wondering why 183 as opposed to 178 enrolled?

finalsubtype_pos <- finalsubtype %>%
  filter(type.inf != "Negative")
# 149 subjects that were positive were represented here in the finalsubtype df
# Which ones of these were ignored in the final set of 142 subject IDs?
# Note that this finalsubtype df is later updated in this script and this may address this question. So I'll hold off on answering this now and see what the updated finalsubtype df looks like.

# Write out this finalsubtype
saveRDS(finalsubtype, file = "Curated Data/Cleaned Data/EMIT_subtypes.RDS")

#### READ in and work with "2016.06.17 1st visit NP swab subtyping II.csv" ####

part2 <- read.csv('UMD_Raw_Data/PCR Data/PCR results/2016.06.17 1st visit NP swab subtyping II.csv', as.is = T)
part2 <- part2 %>% 
  select(-Well, -Well.Type, -Threshold..dR.)

# Number of rows in part2 file
print(nrow(part2))
# Number of columns in part2 file
print(ncol(part2))

# Remove the rows within the data (e.g. reference, calibration, etc.) not needed for analysis
n2 <- part2 %>% 
  filter(Ct..dR. != 'Reference') %>% 
  filter(!grepl('[A-Z]_[A-Z]', Well.Name)) %>% 
  filter(!grepl('PosControl', Well.Name)) 

# Generate a new column for subject_ID which can be obtained from sample_ID
n2 <- n2 %>% 
  mutate(subject_id = gsub('_[1-9]*_[0-Z]*[0-Z]*[0-Z]*[0-Z]*', '', Well.Name)) %>% 
  mutate(subject_id = gsub('nf[0-Z]*', '', subject_id)) %>% 
  mutate(subject_id = gsub('dm[0-Z]*', '', subject_id)) %>%
  mutate(types = gsub('[0-9]*_[0-9]*', '', Well.Name))
n2$subject_id <- as.numeric(n2$subject_id)
n2 <- n2 %>% 
  arrange(subject_id)

#### Nothing ever gets done with the above n2 object. Is this correct? ####
# It looks like all of the data from this n2 df gets incorporated elsewhere - probably from other pcr assays.
# The borderline H3 detection (ct >41) for 226_7 (NPS) is not reflected in the final data and instead is marked as Unsubtypable A

#### Merge the PCR data with sample virus subtype ####

flu.types <- readRDS("Curated Data/Cleaned Data/EMIT_subtypes.RDS")
flu.types <- select(flu.types, subject.id, type.inf)
# Remember here from above that this df has data from 183 subjects, 149 of whom are positive. 

# Include subtype (from the flu.types df, formerly an object called finalsubtype) with the total.pcr file
# Note that the flu.types file has 156 subject IDs (presumably the enrolled and flu positive individuals)
# But how many subject IDs does the total.pcr file have?
total.pcr_subjectID_count <- total.pcr %>%
  group_by(subject.id) %>%
  summarise(count = n())
# Note that there are actually 202 subjectIDs here in this total.pcr file!

# Were 46 of these 202 with pcr data negative and that's how we got the 156 positive cases? 
# There are only 183 subjects included in the subtype analysis, 149 of which are positive.
# Do the 149 positives match the 156 positive cases that we ultimately get?
# If so, then what happened to the other 7 that are positive but not included in the subtype df?

includetype <- inner_join(total.pcr, flu.types, by = "subject.id")
includetype <- rename(includetype, sample.id = Well.Name)

# The includetype df has 1713 observations
includetype_subjectID_check <- includetype %>%
  group_by(subject.id) %>%
  filter(!is.na(subject.id)) %>%
  count()
# Those 1713 observations are made of data from 178 individuals (the 178 enrolled)

includetype1 <- includetype %>% 
  filter(type.inf == 'Negative')
# The includetype1 df has 407 observations

# However, we can already look at this set and see that what might have been classified previously as negative sometimes has a positive for an aerosol sample, or on a subsequent day of testing. We see lots of observations in this set of 407 that have ct values!
# We need to check to see what happens to these data moving forward. 
# The "negative" df generated below addresses this!

# How many subject IDs are in this includetype1 df?
includetype1_subjectID_check <- includetype1 %>%
  group_by(subject.id) %>%
  filter(!is.na(subject.id)) %>%
  count()
print(nrow(includetype1_subjectID_check))
# Those 407 observations are made of data from 29 individuals
# It looks like many of these "subtype negative" but positive pcr individuals will later get their subtype status updated based on the pcr data.

#### Merge PCR data with sample type ####

allsamples <- readRDS("Curated Data/Cleaned Data/EMIT_samples.cc.RDS")
sampletype <- allsamples %>% 
  select(subject.id, sample.id, sample.type)
# This includes 1938 observations, but from how many subject IDs? Let's check.
sampletype_subjectID_check <- sampletype %>%
  group_by(subject.id) %>%
  filter(!is.na(subject.id)) %>%
  count()
print(nrow(sampletype_subjectID_check))
# There were 355 subject IDs here that produce the 1938 observations. 

# It looks like there are some samples where the sample type is missing - need to fill these in!
sampletype_na_check <- allsamples %>% 
  select(subject.id, sample.id, sample.type) %>%
  filter(is.na(sample.type))
# This shows that there are 19 instances (samples) where we are missing the sample type
# This could have implications for downstream data cleaning steps when we merge with this sampletype df!

# Find out the negative virus subtype but have positive pcr results. From these we can define the virus subtypes of some of the negative cases
# I don't understand what the above comment is trying to get at. If a case was negative, then there would be no subtype, correct?
# Unless we are actually talking about negative aerosol samples, which seems to make more sense, since the includetype df includes pcr data only for the aerosol samples
negative <- includetype %>% 
  filter(type.inf == 'Negative' & Ct..dRn. > 0) %>% 
  select(Experiment, sample.id, Ct..dRn., subject.id, type, type.inf) %>%
  inner_join(sampletype, by = c("subject.id", "sample.id"))
print(nrow(negative))
# This gives 27 observations from 9 subject IDs (hmm, really?)
# In this instance the inner_join with the sampletype df does not change the resulting "negative" df - still 27 observations from 9 subject IDs

# Subjects with negative sample type but positive pcr results from either GII sample or 2rd/3rd NP swab.
print(negative)

saveRDS(negative, "Curated Data/Cleaned Data/negative subtype sample with positive pcr.RDS")

## Based on the pcr results from GII samples or 2rd/3rd np, we have modified a few subjects' subtype

#### READ in "negative subtype sample with positive pcr.RDS" ####

updatetype <- readRDS("Curated Data/Cleaned Data/negative subtype sample with positive pcr.RDS")

updatetype1 <- updatetype %>% 
  select(subject.id, type) %>% 
  distinct(subject.id, type, .keep_all = TRUE)

updatetype2 <- updatetype1 %>% 
  filter(type == 'A')

updatetype3 <- updatetype1 %>% 
  filter(type == 'B')

finalsubtype$type.inf[finalsubtype$subject.id == 52] <- 'B'
finalsubtype$type.inf[finalsubtype$subject.id == 58] <- 'B'
finalsubtype$type.inf[finalsubtype$subject.id == 105] <- 'Unsubtypable A' 
finalsubtype$type.inf[finalsubtype$subject.id == 223] <- 'B and unsubtypable A'
finalsubtype$type.inf[finalsubtype$subject.id == 226] <- 'Unsubtypable A' 
finalsubtype$type.inf[finalsubtype$subject.id == 327] <- 'B'
finalsubtype$type.inf[finalsubtype$subject.id == 231] <- 'B'
finalsubtype$type.inf[finalsubtype$subject.id == 329] <- 'B'
finalsubtype$type.inf[finalsubtype$subject.id == 365] <- 'B'

# Note: preveiously Jing had commented out the lines that gave subject IDs 52 and 58 the subtype B designation
# However, there is no evidence to do this, and the pcr files do show that these subject had samples that were positive for flu B
# Jing also had indicated that 223 was a B type of influenza, however, this is not true - 223 is actually type A and B
# Perhaps this was an error of misreading the updatetype df?

#### READ in "EMIT_samples.cc.RDS" ####

enrollcheck <- readRDS('Curated Data/Cleaned Data/EMIT_samples.cc.RDS') %>% 
  select(subject.id, enrolled) %>% 
  distinct(subject.id, enrolled, .keep_all = TRUE) %>% 
  filter(enrolled == TRUE)

finalenrolltype <- finalsubtype %>%
  semi_join(enrollcheck, by = 'subject.id') %>%
  arrange(subject.id) %>%
  select(subject.id, type.inf)

saveRDS(finalenrolltype, "Curated Data/Cleaned Data/EMIT_subtypes_enrolled.RDS")

finalenrolledpositive <- finalenrolltype %>% 
  filter(!type.inf == 'Negative')
# This shows that there were actually 158 positive cases (when we added the 2 flu B's on that were previously commented out: SID 52 and 58)

saveRDS(finalenrolledpositive, "Curated Data/Cleaned Data/EMIT_subtypes_enrolled_positive.RDS")

negative <- finalenrolltype %>% 
  filter(type.inf == 'Negative')

h3n2 <- finalenrolltype %>% 
  filter(type.inf == 'H3N2')

B <- finalenrolltype %>% 
  filter(type.inf == 'B')

Pandemic.H1 <- finalenrolltype %>% 
  filter(type.inf == 'Pandemic H1')

#### Check that all the negative subjects do not have any positive GII or 2nd/3rd np positive PCR samples ####

flu.typesenroll <- readRDS("Curated Data/Cleaned Data/EMIT_subtypes_enrolled.RDS")
flu.typesenroll <- flu.typesenroll %>%
  select(subject.id, type.inf)

# Add the type of infection to the total.pcr df and create new df called includetypeenroll
includetypeenroll <- total.pcr %>%
  inner_join(flu.typesenroll, by = "subject.id") %>%
  rename(sample.id = Well.Name)
# This gives a df with 1,713 observations that should be from the 178 enrolled participants

# Let's check to make sure that these 1,713 observations do indeed come from 178 enrolled participants
includetypeenroll_subjectID_check <- includetypeenroll %>%
  group_by(subject.id) %>%
  filter(!is.na(subject.id)) %>%
  count()
# This is correct! We see that there is data on 178 enrolled participants

includetypeenrollneg <- includetypeenroll %>% 
  filter(type.inf == 'Negative') %>%
  inner_join(allsamples, by = c("subject.id", "sample.id")) %>% 
  filter(focus.ct > 0)
# 50, 234, 306 are negative cases with all samples negative for PCR, but 50_3, 234_3, 306_3 are positive for focus assay

#### READ in "EMIT_subtypes_enrolled_positive.RDS" ####

flu.typepositive <- readRDS("Curated Data/Cleaned Data/EMIT_subtypes_enrolled_positive.RDS")

#### Check the unmatched cases between sample virus subtype and PCR results ####

# Add the infection types from the positive cases to the total.pcr df and add the sample types to each of the sample ids
includetypepostive <- total.pcr %>%
  inner_join(flu.typepositive, by = "subject.id") %>% 
  rename(sample.id = Well.Name) %>%
  inner_join(sampletype, by = c("subject.id", "sample.id"))
# *** Note: using an inner_join here for the sampletype eliminates 2 observations in the includetypepositive df that would otherwise be there were the sampletype df complete and did not contain NAs - there are 19 observations in the sampletype df that are missing their sample type (i.e., NPS, Fine, Coarse, Throat, Ant Nares, etc.)
# I may be able to locate these samples in the freezer and determine the sample types. 

# If the inner_join with the sampletype df with missing sample types is used, then we get a "includetypepositive" df with 1461 observations
# How many subject IDs contribute to these 1461 observations.
# *Remember that if we fix the sampletype file (i.e., fill in the 19 missing sample types), there is potential for 1463 observations
includetypepositive_subjectID_check <- includetypepostive %>%
  group_by(subject.id) %>%
  filter(!is.na(subject.id)) %>%
  count()
print(nrow(includetypepositive_subjectID_check))
# There are 158 subjectIDs represented in this includetypepostive df

# Let's check to see if there are any intances where the type is B but the subtype (type.inf variable) is not B

unmatched1 <- includetypepostive %>% 
  filter(type == 'B' & !grepl('B', type.inf))
print(unmatched1)
# Even though there are a bunch of observations here that are labeled type == 'B' but type.inf not 'B', none of these observations have pcr evidence for B virus.
# Thus, the type.inf appears to be not incorrect for each of these observations. We would have to check to see that there was pcr evidence for the type of infection labeled by type.inf in order to prove that these type.inf labels were correct. 

# Let's also check to see if there are any instances where the type is A but the subtype (type.inf variable) is not A 
unmatched2 <- includetypepostive %>%
  filter( type == 'A' &  type.inf == 'B')

# Similarly, there were no instances where there was a type.inf of B where any samples were pcr positive for A virus. 

# Does this mean we can eliminate these observations in the unmatched1 and unmatched2 dfs? Yes.

includetypepostiveupdate <- includetypepostive %>%
  anti_join(unmatched1, by = c("Experiment", "sample.id", "Ct..dRn.", "copies.in", "subject.id", "type", "copy.num", "date", "cfactor", "virus.copies", "type.inf", "sample.type")) %>% 
  anti_join(unmatched2 , by = c("Experiment", "sample.id", "Ct..dRn.", "copies.in", "subject.id", "type", "copy.num", "date", "cfactor", "virus.copies", "type.inf", "sample.type"))

# SubjectID check on includetypepositiveupdate

includetypepostiveupdate_subjectID_check <- includetypepostiveupdate %>%
  group_by(subject.id) %>%
  filter(!is.na(subject.id)) %>%
  count()
# We see that we have observations for all 158 enrolled and positive subjects
# We notice that subjectID 52 has extraneous observations because there are 2 identical pcr experiments that have slightly different names
# One is called: "15. 2015.07.29.013.243 2012-13 samples FluA_B and FluA PCR"
# Another is called: "17. 2015.7.29.013.243 2012-13 samples FluA_B and FluA PCR"

# *** We can go back and fix this in the pcr files

#### Pick pos fluA, pos fluB, join, assign copy number ####

# Pick out all the PCR with A assay results
allpcrA <- includetypepostiveupdate %>% 
  filter(type == 'A') %>%
  rename(virus.copiesA = virus.copies, typeA = type) %>%
  rename(CtA = Ct..dRn.)

# Pick out all the PCR with B assay results
allpcrB <- includetypepostiveupdate %>% 
  filter(type == 'B') %>%
  rename(virus.copiesB = virus.copies, typeB = type) %>%
  rename(CtB = Ct..dRn.)

# Join both A and B assay, and also add sample type in the data list, assign the final RNA copies number for each sample type
allPCR <- allpcrA %>%
  full_join(allpcrB, by = c('subject.id', 'Experiment', 'type.inf', 'sample.id', 'sample.type', "copies.in", 'copy.num', 'date', 'cfactor')) %>%
  arrange(subject.id) %>% 
  filter(!sample.type == 'Throat Swab')

## DATA EDITING (dilution factors different for a few samples) ##
# Flagged samples all NP swabs, dilution factor is 50
# 66_7 120_7 184_8 188_7 189_7 192_7 196_7 262_7 277_7 284_12 284_7 296_12 296_7

allPCR1 <- allPCR %>% 
  filter(sample.id == '66_7' | sample.id == '120_7' | sample.id == '184_8' | sample.id == '188_7' | sample.id == '189_7' | sample.id == '192_7' | 
           sample.id == '196_7' | sample.id == '262_7' | sample.id == '277_7' | sample.id == '284_12' | sample.id == '284_7' | 
           sample.id == '296_12'| sample.id=='296_7')

allPCR2 <- allPCR %>%
  anti_join(allPCR1)

allPCR3 <- allPCR1 %>% 
  mutate(final.copiesA = virus.copiesA*50, final.copiesB = virus.copiesB*50)

allPCR4 <- allPCR2 %>% 
  filter(!sample.type == 'Nasopharyngeal swab') %>% 
  mutate(final.copiesA = virus.copiesA*25, final.copiesB = virus.copiesB*25)

allPCR5 <- allPCR2 %>% 
  filter(sample.type == 'Nasopharyngeal swab') %>%
  mutate(final.copiesA = virus.copiesA*100, final.copiesB = virus.copiesB*100)

allPCRtotal <- rbind(allPCR3, allPCR4, allPCR5) %>% 
  ungroup

allPCR.A <- allPCRtotal %>% 
  filter(typeA == 'A') %>%
  select(-CtB, -virus.copiesB, -typeB, -final.copiesB) %>%
  rename(Ct = CtA, virus.copies = virus.copiesA, type = typeA, final.copies = final.copiesA)

allPCR.B <- allPCRtotal %>% 
  filter(typeB == 'B') %>%
  select(-CtA, -virus.copiesA, -typeA, -final.copiesA) %>%
  rename(Ct = CtB, virus.copies = virus.copiesB, type = typeB, final.copies = final.copiesB)

# 2. merge the seperated FILE FOR A and B back together, successfully make,- one row for one subject.id
allPCRfinal <- rbind(allPCR.A, allPCR.B) %>% 
  ungroup

allPCRfinal_subjectID_check <- allPCRfinal %>%
  group_by(subject.id) %>%
  count()
print(nrow(allPCRfinal_subjectID_check))
# Data exists here for the 158 positive enrolled cases

saveRDS(allPCRfinal, "Curated Data/Cleaned Data/pcr data for gii and 2nd&3rd np.RDS")

#### Check to make sure no repeats ####

no.repeats.pos <- allPCRfinal %>% 
  group_by(sample.id) %>%
  summarise(n = n())
# This is useful, there should be 2 repeats per sample (becasue samples were run in duplicate)
# However, in cases where a sample was positive for both fluA and B, then there will be 4 pcr data repeats per sample (a duplicate for flu A assay and a duplicate for flu B assay)
# However, in the case where we ran out of sample and went to a new tube and the new tube had a slighly different label (for example, samples with subscripts _1 and _7 are often NPS for the first day of sampling for a particular subject)

# Subjects with known sample type should only have 2 PCR result
no.repeats.pos.1 <- no.repeats.pos %>% 
  filter(n == 1) %>% 
  left_join(allPCR, by = 'sample.id') %>% 
  select(Experiment, sample.id, CtA, copies.in, typeA, type.inf, CtB, typeB)
# Subjects with known sample type have only one PCR result (Hmm - I'm don't understand this)
print(no.repeats.pos.1)

no.repeats.pos.2 <- no.repeats.pos %>% 
  filter(n == 4) %>%
  left_join(allPCR, by = 'sample.id') %>%
  select(Experiment, sample.id, CtA, copies.in, typeA, type.inf, CtB, typeB)
# Subjects with known sample type have two PCR result
print(no.repeats.pos.2)

## This is interesting because it shows when a sample was repeated on a pcr assay (not talking about duplicates here, but rather repeat assays). 
# Is there a rule about which of the repeats to keep and which to exclude? For example, should the most recent assay be taken and the other excluded?
# Does this ever get addressed later in any code?

#### **** Using Script: Jing and Dr. Milton's "1st np swab quantity.R" **** ####

###
## Original file information:

# Author: Jing Yan & Don Milton
# Date: September 17, 2015
# Revision Date: Jun 30, 2016
# Title: 1st np swab quantity.R
# Purpose: To sort the data files from the lab (Michael Grantham) PCR data for 1st np samples 

# Input files:
# InputFiles_UMD/PCR results/2016.06.17 1st visit NP swab FluA quant.csv
# InputFiles_UMD/PCR results/2016.06.17 1st visit NP swab FluB quant.csv
# Output files: R_output

### Procedures:

# 1. Sort out first NP swab PCR results(includes A and B, combine the two parts after sorting)
# 2. For the first NP swab, seperate the subjects with flu A infection, flu B infection, 
#    negative on the swab and dual infection of A and B
# 3. Question: How to treat the samples with multiple PCR results? 

## Method 1
# 1 obs---use
# >=2 obs---take mean(if one is No Ct, Tobit fitted value)

## Method 2
# Take a tobit model(obs~sample_id)---get fitted data for all the subjects
# (check if fitted value match with the method1 values)

## Method 3 
# Tobit(obs~sample.type+subject.id)---get fitted data for all the subjects (may be different from Method 1)

# Question about the above notes (the notes above are from Jing) = how is Method 2 different from Method 3? They both say that they are meant to "get fitted data for all the subjects". Is Method 3 trying to get fitted data for all of the samples?

###
#### READ in and work with "EMIT_subtypes_enrolled_positive.RDS" ####
## Note that this df is the product of other script and in Jing's original setup, was saved to R_output
# The original script from Jing that produced this df was: "subtype analysis.R"
## However, in this new setup, we have saved the EMIT_subtypes_enrolled_positive.RDS in: 
# ... EMIT/EMIT_Data_Analysis_Jake/EMIT_UMD_Natural_Infection/Curated Data/Cleaned Data

flu.types <- readRDS("Curated Data/Cleaned Data/EMIT_subtypes_enrolled_positive.RDS")
flu.types <- select(flu.types, subject.id, type.inf)

#### READ in "fluA_1np_calibration.RDS" ####
npAcali <- readRDS("Curated Data/Cleaned Data/fluA_1np_calibration.RDS")

#### READ in "2016.08.05 1st visit NP swab FluA quant.csv" ####
npA <- read.csv('UMD_Raw_Data/PCR Data/PCR results/2016.08.05 1st visit NP swab FluA quant.csv', as.is = T)

#### Work with the above 2 dfs ####
names(npA)
npA <- npA %>%
  filter(Ct..dRn. != 'Reference') %>%
  filter(!grepl('_Low', Well.Name)) %>%
  filter(!grepl('NTC', Well.Name)) %>%
  filter(!grepl('Standard', Well.Type)) %>%
  filter(!grepl('High', Well.Name)) %>%
  filter(grepl('_', Well.Name)) %>%
  select(-Well, -Well.Type, -Threshold..dRn.) %>%
  rename(copies.in = Quantity..copies.)

# Number of rows in part1 first NP influenza A PCR data
nrow(npA)
# Number of column in part1 first NP influenza A PCR data
ncol(npA)

npA <- npA %>% 
  mutate(subject.id = gsub('_[1-9]*_[0-Z]*', '', Well.Name)) %>% 
  mutate(subject.id = gsub('_A', '', subject.id)) %>% 
  mutate(subject.id = gsub('nfA', '', subject.id))
npA[npA == "No Ct"] <- ''
npA$copies.in <- as.numeric(npA$copies.in)
npA <- npA %>%
  mutate(copy.num = copies.in*100*250) %>%
  arrange(subject.id)
npA$result.type <- 'A'
npA <- npA %>%
  mutate(sample.id = gsub('_A', '', Well.Name)) %>% 
  mutate(sample.id = gsub('_InfA', '', sample.id))
npA$sample.id[npA$sample.id == '118_1_1'] <- '118_1' 
npA$sample.id[npA$sample.id == '69_1_1'] <- '69_1' 
npA$sample.id[npA$sample.id == '70_1_1'] <- '70_1' 
npA$sample.id[npA$sample.id == '210_1_1'] <- '210_1'
npA$sample.id[npA$sample.id == '339_1_1'] <- '339_1'
npA$sample.id[npA$sample.id == '356_1_1'] <- '356_1'
npA <- npA %>%
  distinct(subject.id, Ct..dRn., .keep_all = TRUE)

# Number of rows in first NP influenza A PCR data
nrow(npA)

# Number of columns in first NP influenza A PCR data
ncol(npA)

npA <- npA %>%
  rename(Ct = Ct..dRn., type = result.type) %>%
  select(-Well.Name) %>%
  mutate(date = gsub('^[0-9]*.', '', Experiment))

npAfinal <- npA %>%
  left_join(npAcali, by = 'date')
# The following samples were done without inter-run calibrators in the experiments, so the adjustment calibrators are missing
npAfinal$cfactor[npAfinal$sample.id == '176_6'] <- 1
npAfinal$cfactor[npAfinal$sample.id == '64_6'] <- 1
npAfinal$virus.copies <- npAfinal$copy.num*npAfinal$cfactor

#### Now working with the flu B calibrations ####

#### READ in "fluB_1np_calibration.RDS" ####
npBcali <- readRDS("Curated Data/Cleaned Data/fluB_1np_calibration.RDS")

#### READ in "2016.08.05 1st visit NP swab FluB quant.csv" ####
npB <- read.csv('UMD_Raw_Data/PCR Data/PCR results/2016.08.05 1st visit NP swab FluB quant.csv', as.is = T)

names(npB)

npB <- npB %>%
  filter(Ct..dRn. != 'Reference') %>%
  filter(!grepl('_Low', Well.Name)) %>%
  filter(!grepl('NTC', Well.Name)) %>%
  filter(!grepl('Standard', Well.Type)) %>%
  filter(!grepl('High', Well.Name)) %>%
  filter(grepl('_', Well.Name)) %>%
  select(-Well, -Well.Type, -Threshold..dRn.) %>% 
  rename(copies.in = Quantity..copies.)

# Number of rows in part1 first NP influenza A PCR data
nrow(npB)
# Number of column in part1 first NP influenza A PCR data
ncol(npB)

npB <- npB %>%
  mutate(subject.id = gsub('_[1-9]*_[0-Z]*', '', Well.Name)) %>%
  mutate(subject.id = gsub('_B', '', subject.id)) %>%
  mutate(subject.id = gsub('nfB', '', subject.id))
npB[npB == "No Ct"] <- ''
npB$copies.in <- as.numeric(npB$copies.in)
npB <- npB %>% 
  mutate(copy.num = copies.in*100*272) %>%
  arrange(subject.id)
npB$result.type <-'B'
npB <- npB %>%
  mutate(sample.id = gsub('_B', '', Well.Name)) %>%
  mutate(sample.id = gsub('_InfB', '', sample.id))
npB$sample.id[npB$sample.id == '118_1_1'] <- '118_1' 
npB$sample.id[npB$sample.id == '210_1_1'] <- '210_1' 
npB$sample.id[npB$sample.id == '339_1_1'] <- '339_1' 
npB$sample.id[npB$sample.id == '356_1_1'] <- '356_1' 
npB$sample.id[npB$sample.id == '69_1_1'] <- '69_1' 
npB$sample.id[npB$sample.id == '70_1_1'] <- '70_1' 
npB <- npB %>%
  distinct(subject.id, Ct..dRn., .keep_all = TRUE)

# Number of rows in first NP influenza A PCR data
nrow(npB)
# Number of columns in first NP influenza A PCR data
ncol(npB)

npB <- npB %>%
  rename(Ct = Ct..dRn., type = result.type) %>%
  select(-Well.Name) %>%
  mutate(date = gsub('^[0-9]*.', '',Experiment))

npBfinal <- npB %>%
  left_join(npBcali, by = 'date')
npBfinal$virus.copies <- npBfinal$copy.num*npBfinal$cfactor

#### Merging together the fluA and fluB data ####

npfirst <- rbind(npAfinal, npBfinal)
npfirst$subject.id = as.integer(npfirst$subject.id)
npfirst <- npfirst[order(npfirst$subject.id), ]

# Check the number of subjectIDs, sampleIDs, and type of infection
npfirst_sampleID_check <- npfirst %>%
  group_by(sample.id, type) %>%
  count()

# npfirst has 378 observations and now let's see how many subject IDs there are

npfirst_subjectID_check <- npfirst %>%
  group_by(subject.id) %>%
  count()
# There are 183 subject IDs here. Most have 2 pcr results, some have more

# Let's look at the subjectIDs with more than 2 pcr results
npfirst_subjectID_check_observations <- npfirst %>%
  group_by(subject.id) %>%
  count() %>%
  filter(n >2)
# There are 11 subjects with multiple day 1 NP swabs that were run for either flu A or flu B. 
# I'm not talking about duplicates in the pcr assay - these are completely new experiments with new dates.
# We will use all of the data, including that from repeat experiments

npfirstpositive <- npfirst %>% 
  inner_join(flu.types, by = c("subject.id"))

npfirstpositive_subjectID_check <- npfirstpositive %>%
  group_by(subject.id) %>%
  count()

# Note: It appears that many of the 183 in the "npfirst" df were actually positive for flu virus, yet these were not included in the "npfirstpositive" df that includes data from the 158 enrolled and positive flu cases.
# We need to do a review to better understand why there were some of these samples that were not included in the final dataset. 

unmatchedfirstnp1 <- npfirstpositive %>% 
  filter( type == 'B' & type.inf == 'H3N2')

unmatchedfirstnp2 <- npfirstpositive %>% 
  filter( type == 'B' & type.inf == 'Pandemic H1')

unmatchedfirstnp3 <- npfirstpositive %>% 
  filter( type == 'B' & type.inf == 'Unsubtypable A')

unmatchedfirstnp4 <- npfirstpositive %>% 
  filter( type == 'B' & type.inf == 'H3N2 and PH1')

unmatchedfirstnp5 <- npfirstpositive %>% 
  filter( type == 'A' &  type.inf == 'B')

unmatchedtotal <- rbind(unmatchedfirstnp1, unmatchedfirstnp2, unmatchedfirstnp3, unmatchedfirstnp4, unmatchedfirstnp5)

# npfirstpositiveupdate <- setdiff(npfirstpositive, unmatchedtotal)
# I actually think what Jing means here is to take the anti_join, to remove the unmatchedtotal observations from the npfirstpositive df. Let's try this and comment out the setdiff command

npfirstpositiveupdate <- npfirstpositive %>%
  anti_join(unmatchedtotal) %>%
  arrange(subject.id)
npfirstpositiveupdate$sample.type <- 'Nasopharyngeal swab'

# Let's look at the number of subject IDs in this df and also check to if any samples have repeat assays
npfirstpositiveupdate_subjectID_check <- npfirstpositiveupdate %>%
  group_by(subject.id) %>%
  summarise(n = n()) %>%
  arrange(n)
# Samples that are duplicated (more than 1 assay) in the "npfirstpositiveupdate" df at least once are: 
# 64, 81, 97, 100, 104, 118, 123, 174, 176, 223, 230, 357, 95
# Need to decide how to deal with these - average?
# However some are negative on one run and positive on separate run - take the positive? the negative?
# Could treat it like when there is one/two positive detections for a replicate on the PCR plate?

# Not sure the following bit of code makes changes the df at all and not really sure of the intended purpose of this bit of code. I will comment it out for now. 
# I actually think the below code is responsible for misclassification of flu virus type (A or B) for samples that had both A and B infection!
# npfirstpositiveupdate <- npfirstpositiveupdate %>%
#   mutate(AorB = ifelse(type.inf == "B", "B",
#                        ifelse(type.inf == "B and unsubtypable A", "A and B",
#                               ifelse(type.inf == "H3N2 and B", "A and B", "A")))) %>%
#   filter(AorB == type | AorB == "A and B") %>%
#   select(-AorB)

saveRDS(npfirstpositiveupdate, "Curated Data/Cleaned Data/EMIT_np_quantity.RDS")

# Seems like there should be 2 PCR results for each of these samples (replicates in the PCR assay) but most of these only have a single pcr copy number. Is this intentional? 


#### Working with the "EMIT_UMD_Final_Data_Summary.R" script ####
# Title: EMIT_UMD_Final_Data_Summary.R
# Moving this file, originally produced by Jing and Don Milton, to git lab
# Also organizing it and cleaning it up some. 

## Original File information

# Author: Jing Yan
# Date: May 03, 2016
# Revision Date: 2016
# Title: final UMD data summary.R
# Purpose: clean data and set files for tobit models
# Input files:pcr data for gii and 2nd&3rd np.RDS
#             EMIT_np_quantity.RDS
#            EMIT_subtypes_enrolled_positive.RDS
#             EMIT_samples.cc.RDS
#          ???EMITClinicalUMD2013.csv
# Output file: tables.txt

#_____________
# Setup all I/O
# setwd('C:/Users/Jing/Box Sync/Box Sync/EMIT/EMIT_Data_Analysis')
# setwd('/Users/dmilton/Box Sync/0_DKM/Lab/Biodefense/EMIT/EMIT_Data_Analysis')
# setwd('/Volumes/Internal RAID Set 1/Box Sync/0_DKM/Lab/Biodefense/EMIT/EMIT_Data_Analysis')
# sink(file="R_output/tables.txt",split=TRUE)
# Out.dir <- "R_output/"
#_____________

pcr1 <- readRDS("Curated Data/Cleaned Data/pcr data for gii and 2nd&3rd np.RDS")

pcr1 <- pcr1 %>% 
  select(-Ct, -copies.in, -copy.num, -date, -cfactor, -virus.copies)

nppcr1 <- readRDS("Curated Data/Cleaned Data/EMIT_np_quantity.RDS")

nppcr1 <- nppcr1 %>% 
  select(-Ct, -copies.in, -copy.num, -date, -cfactor) %>% 
  rename(final.copies = virus.copies) # This seems to be incorrect. To get from virus.copies to final.copies one must multiply virus.copies by the dilution factor. For NP swabs the dilution factor was 100 (exept for a few NPS where 100ul instead of 50 ul was used to extract, in which case the dilution factor was 50).

# Merge the preceding 2 dfs
allpcr <- rbind(pcr1, nppcr1) %>% 
  ungroup

# Check number of subjectIDs here
allpcr_subjectID_check <- allpcr %>%
  group_by(subject.id) %>%
  count()
print(nrow(allpcr_subjectID_check))

# 114_1, 127_1, 335_1 (These three subjects were enrolled on the second visit, these three samples are from their first screen visit) 
# I presume we are removing these because they don't have aerosol samples on their first day of sample collection, while they do have NPS data?
# Subject 333 was removed, because the coarse sample was lost
allpcr <- allpcr[!(allpcr$sample.id == "114_1"| allpcr$sample.id == "127_1" | allpcr$sample.id == "335_1" | allpcr$sample.id == "333_3" | allpcr$sample.id == "333_1"), ]

# Now that some data has been removed, again check number of subjectIDs here
allpcr_subjectID_check <- allpcr %>%
  group_by(subject.id) %>%
  count()
print(nrow(allpcr_subjectID_check))

# I'm confused why the data for 114, 127, and 335 was simply removed for the first visit, but not the subsequent visits - wasn't the point of removing this data to eliminate the full pcr data for these subjects, because they don't have aerosol data on the first day of study where they do have NPS data? This is because they had data for dpo = 0 so that dpo = 0 visit was removed, while other data that was for dpo=1 or dpo=2 or dpo=3 was kept.

## Read EMIT_samples.cc.RDS
all <- readRDS("Curated Data/Cleaned Data/EMIT_samples.cc.RDS")

subgroup <- all %>% 
  select(subject.id, sample.id, date.visit, date_on_sx)
subgroup$date.visit <- as.Date(subgroup$date.visit)

passagefocus <- all %>% 
  select(subject.id, sample.id, date.visit, passpos, validp, focus.ct)

coughgiivisits <- all %>% 
  select(subject.id, sample.id, date.visit, cough_number, sneeze_number, sample.type, g2.run)
# Does this cough_number have all the data in it?
# There seems to be quite a few NA's in there!
# May need to find a cough_number variable that has all of the data in it.
# Note: that in earlier dfs that were used to create this object, there were missing cough data with notes that the audio recordings should be used to identify cough number for quite a few GII sampling instances. 

# Merge and clean
pcrgiivisit <- allpcr %>%
  inner_join(coughgiivisits, by = c('subject.id', 'sample.id', 'sample.type'))
pcrgiivisit$date.visit <- as.Date(pcrgiivisit$date.visit, "%m/%d/%Y")

pcrgiivisit_subjectID_check <- pcrgiivisit %>%
  group_by(subject.id) %>%
  count()

## Read daypostonset.csv 
# (note: it is unclear how this file was created. It was found in Jing's R_output folder so presumably it was created with some script - unfortunately a search of files yielded none that might have produced this file)
daypick1 <- read.csv("UMD_Raw_Data/REDCAP/daypostonset.csv")

# New plan (Januaray 29, 2019) regarding use of daypostonset data
# It looks like the daypostonset.csv file that was in the raw data only has observations from 149 subject IDs and this is actually a major source of lost subject ID data between the 158 enrolled and positive cases and the 142 that made it into Jing's PNAS_df. 

# Let's go back to the samples.cc df to recreate a new dayposonset df that has the data for all 178 screened individuals.

samples.cc_dpo <- samples.cc %>%
  distinct(subject.id, date.visit, date_on_sx)
samples.cc_dpo$date.visit <- as.Date(samples.cc_dpo$date.visit)
samples.cc_dpo$date_on_sx <- as.Date(samples.cc_dpo$date_on_sx)
samples.cc_dpo <- samples.cc_dpo %>%
  mutate(dpo = date.visit - date_on_sx) %>%
  arrange(subject.id)

# How many subject IDs is this for?
samples.cc_dpo_subjectID_check <- samples.cc_dpo %>%
  group_by(subject.id) %>%
  count()

# Need to add in the dpo where they are NA based on the date of symptom onset!

sub <- unique(samples.cc_dpo$subject.id)
samples.cc_dpo_full <- samples.cc_dpo %>%
  filter(is.na(subject.id))

for (i in 1:length(sub)) {
  subid <- sub[i]
  temp <- samples.cc_dpo[samples.cc_dpo$subject.id == subid, ]
  for (j in 1:(nrow(temp))) {
    temp$dpo[j] = temp$date.visit[j] - temp$date_on_sx[1]
    samples.cc_dpo_full <- rbind(samples.cc_dpo_full, temp)
  }
}

# This loop works ok, but there is a glitch with the NAs
# To overcome this, we will filter out where dpo = NA and then take only the unique rows (to eliminate duplication that occurred in the loop)
samples.cc_dpo <- samples.cc_dpo_full %>%
  filter(!is.na(dpo)) %>%
  distinct(subject.id, dpo, .keep_all = TRUE)

# Note there were a few times where there was no date_on_sx and so dpo could not be determined.
# This set of individuals was excluded from this version of samples.cc_dpo

# check the subjectIDs remaining in samples.cc_dpo
samples.cc_dpo_subjectID_check <- samples.cc_dpo %>%
  group_by(subject.id) %>%
  count()
# There are 331 subjectIDs here - this means that we lost 24 that didn't have a date_on_sx recorded
# Hopefully these 24 subjectIDs are part of the unenrolled group, however we will be able to check this.

# Make the samples.cc_dpo df take on the daypick1 df
daypick1 <- samples.cc_dpo

## old version
# daypick <- daypick1 %>%
#   distinct(subject.id, date.visit, dpo, .keep_all = TRUE) %>%
#   arrange(subject.id, dpo)
# 
# daypickfinal <- daypick %>% 
#   filter(!(dpo == 0 | dpo == 4 |dpo == 7))
# daypickfinal$date.visit <- as.Date(daypickfinal$date.visit, format = "%m/%d/%Y")
# daypickfinal$date_on_sx <- as.Date(daypickfinal$date_on_sx, format = "%m/%d/%Y")
## end old version

daypickfinal_full <- daypick1
daypickfinal_full$dpo <- as.numeric(daypickfinal_full$dpo)
daypickfinal_full <- daypickfinal_full %>%
  arrange(subject.id, dpo)

# SubjectID count check for daypickfinal_full
daypickfinal_full_subjectID_check <- daypickfinal_full %>%
  group_by(subject.id) %>%
  count()
print(nrow(daypickfinal_full_subjectID_check))
# 331 subjectIDs
# But how many sampling instances?
daypickfinal_full_sampling_instance_check <- daypickfinal_full %>%
  distinct(subject.id, date.visit)
print(nrow(daypickfinal_full_sampling_instance_check))
# 440 sampling instances when we use all of the daypickfinal data (inclusive of dpo day)

# But when we restrict the dpo day to just dpo == 1, 2, or 3, then how many sampling instances do we get?

# Jing's final PNAS dataset requires this daypickfinal to only include the dpo=1-3, however we left the daypickfinal_full version above to help recreate and understand all the places we whittled down subjectIDs and/or sampling instances
daypickfinal_dpo1_through_3 <- daypick1 %>% 
  filter(dpo == 1 | dpo == 2 | dpo == 3)
daypickfinal_dpo1_through_3$dpo <- as.numeric(daypickfinal_dpo1_through_3$dpo)
daypickfinal_dpo1_through_3 <- daypickfinal_dpo1_through_3 %>%
  arrange(subject.id, dpo)

# SubjectID count check for daypickfinal_dpo1_through_3 (before merge)
daypickfinal_dpo1_through_3_subjectID_check <- daypickfinal_dpo1_through_3 %>%
  group_by(subject.id) %>%
  count()
print(nrow(daypickfinal_dpo1_through_3_subjectID_check))

# This new daypickfinal_dpo1_through_3 df has 309 subjectIDs (the filter command to select only those with dpo == 1, 2, or 3 eliminated 22 subjectID entries from daypickfinal_full)
# But how many sampling instances do we get now with the daypickfinal_dpo1_through_3 df?
daypickfinal_dpo1_through_3_sampling_instance_check <- daypickfinal_dpo1_through_3 %>%
  distinct(subject.id, date.visit)
print(nrow(daypickfinal_dpo1_through_3_sampling_instance_check))
# we get 399 sampling instances when restricting the dpo df to dpo=1-3

# Merge and clean
withdpo <- pcrgiivisit %>%
  inner_join(daypickfinal_dpo1_through_3, by = c('subject.id', 'date.visit'))

# Subject ID count check for withdpo
withdpo_subjectID_check <- withdpo %>%
  group_by(subject.id) %>%
  count()
print(nrow(withdpo_subjectID_check))

# This merge reduces the number of subjectIDs from 157 in the pcrgiivisit df to 148 in the withdpo df (the daypickfinal_dpo1_through_3 df had 309)
# Let's examine to see how we lost these 9 subjectIDs
# Any subjectIDs in pcrgiivisit that weren't in daypickfinal_dpo1_through_3?
pcrgiivisit_but_not_in_daypickfinal_dpo1_through_3 <- pcrgiivisit %>%
  anti_join(daypickfinal_dpo1_through_3, by = c('subject.id')) %>%
  group_by(subject.id) %>%
  count()
print(nrow(pcrgiivisit_but_not_in_daypickfinal_dpo1_through_3))
# In checking the dpo df before excluding for dpo == less than 1 or greater than 3, we see that:
# 122 only had a dpo 4 so was excluded
# 166 only had a dpo 4 so was excluded
# 249 only had a dpo 7 so was excluded
# 299 only had a dpo 4 so was excluded
# 302 only had a dpo 4 so was excluded
# 318 only had a dpo 4 so was excluded
# 327 only had a dpo 0 so was excluded
# 350 only had a dpo 0 so was excluded

# But what about the 9th missing subjectID?
# Let's compare the merged df with the pcrgiivisit df to see who is missing
pcrgiivisit_but_not_in_withdpo <- pcrgiivisit %>%
  anti_join(withdpo, by = "subject.id") %>%
  group_by(subject.id) %>%
  count()
print(nrow(pcrgiivisit_but_not_in_withdpo))

ninth_missing_subjectID <- pcrgiivisit_but_not_in_withdpo %>%
  anti_join(pcrgiivisit_but_not_in_daypickfinal_dpo1_through_3)
# This shows that 114 was the other missing subjectID
# Why was 114 missing? Was it missing from the daypickfinal_dpo1_through_3 df?

record_for_114 <- pcrgiivisit %>%
  filter(subject.id == 114)
# We see that 114 had positive pcr detections on dpo 4 but there were no pcr records for 114 from the only other day of sample collection, which was on dpo 1. There shoud be more pcr records for 114 from dpo 1. Even if these were all negative, this case could be included in the final PNAS_df just like 174 was! Otherwise there is inconsistency in the analytical inclusion criteria!

withdpo <- withdpo %>% 
  filter(!(is.na(cough_number)))
# But does this cough_number include all the data from the recordings? We assume so.

# Subject ID count check for withdpo
withdpo_subjectID_check <- withdpo %>%
  group_by(subject.id) %>%
  count()
print(nrow(withdpo_subjectID_check))
# Thus - there were 2 records among the 148 that didn't have cough data so those were removed and now we are down to 146.

## Jing says:
# Remove experiment 9. 2015.06.19.014.74 2012-2013 Samples PCR Flu B 
# Then we removed subjects 322 and 337 and remove 182 second visit
## End of Jing comment

# But, what is the explanation for this? We solved the case - see note below - but it was due to bad interrun calibrator in the pcr assay.

withdpo_322_record <- withdpo %>%
  filter(subject.id == 322)
# Only thing out of place is that the first sample was taken on dpo=0, but this shouldn't be reason to exclude the dpo=1 samples, remaining here, from inclusion!
# I haven't seen any reason to support exclusion of this case. 

withdpo_337_record <- withdpo %>%
  filter(subject.id == 337)
# I haven't seen any reason to support exclusion of this case. 
# Is it because there is only a dpo=3 record here?

withdpo_182_record <- withdpo %>%
  filter(subject.id == 182)
# I haven't seen any reason to support exclusion of dpo = 2 data for this subject. 
# Is it because there is no dpo = 1 data here? If so, then why keep the dpo=3 data?

withdpo <- withdpo[!(withdpo$subject.id == 322 | withdpo$subject.id == 337), ]
withdpo <- withdpo[!(withdpo$subject.id == 182 & withdpo$g2.run == 2), ]

# Subject ID count check for withdpo
withdpo_subjectID_check <- withdpo %>%
  group_by(subject.id) %>%
  count()

# After check some more notes, we see that 322, 337 and the second day of pcr data for 182 were excluded because the interrun-calibrator in the pcr assay was not working properly. 

clinical_in_file <- 'UMD_Raw_Data/REDCAP/EMITClinicalUMD2013.csv'
clinical_umd <- read.csv(clinical_in_file)

# Check out the format of the date_visit variable
clinical_umd_date <- clinical_umd %>%
  select(date_visit)

clinical_umd$date_visit <- as.Date(clinical_umd$date_visit, format = "%m/%d/%y")

clinical_umd1 <- clinical_umd %>% 
  select(field_subj_id, body_temp, date_visit, nose_run, nose_stuf, sneeze, throat_sr, earache, malaise, headache, mj_ache, sw_fever_chill, lymph_node, chest_tight, sob, cough, fluvac_cur, sex, asthma, anitviral_24h, smoke) %>%
  rename(subject.id = field_subj_id, date.visit = date_visit)

clinical_umd_1 <- clinical_umd1 %>%
  select(subject.id, body_temp, date.visit) %>%
  filter(!is.na(body_temp))

# Subject ID count check for "clinical_umd_1" df
clinical_umd_1_subjectID_check <- clinical_umd_1 %>%
  group_by(subject.id) %>%
  count()
print(nrow(clinical_umd_1_subjectID_check))

clinical_umd_2 <- clinical_umd1 %>% 
  select(subject.id, date.visit, nose_run, nose_stuf, sneeze, throat_sr, earache, malaise, headache, mj_ache, sw_fever_chill, lymph_node, chest_tight, sob, cough) %>% 
  filter(!is.na(nose_run))

# Subject ID count check for "clinical_umd_2" df
clinical_umd_2_subjectID_check <- clinical_umd_2 %>%
  group_by(subject.id) %>%
  count()
print(nrow(clinical_umd_2_subjectID_check))

clinical_umd_3 <- clinical_umd1 %>% 
  select(subject.id, asthma) %>% 
  filter(!is.na(asthma))
# Save this object for later merge

# Subject ID count check for "clinical_umd_3" df
clinical_umd_3_subjectID_check <- clinical_umd_3 %>%
  group_by(subject.id) %>%
  count()
print(nrow(clinical_umd_3_subjectID_check))

# Read in new df
subtype <- readRDS("Curated Data/Cleaned Data/EMIT_subtypes_enrolled_positive.RDS")

# Merge and clean
comfirmcases <- clinical_umd_2 %>%
  inner_join(subtype, by = 'subject.id')
# Subject ID count check
comfirmcases_subjectID_check <- comfirmcases %>%
  group_by(subject.id) %>%
  count()

comfirmcases1 <- comfirmcases %>% 
  mutate(upper_sym = nose_run + nose_stuf + sneeze + throat_sr + earache) %>% 
  mutate(lower_sym = chest_tight + sob + cough) %>% 
  mutate(systemic_sym = malaise + headache + mj_ache + lymph_node + sw_fever_chill)

# Merge
comfirmcases2 <- comfirmcases1 %>%
  inner_join(subgroup, by = c('subject.id', 'date.visit'))
# Subject ID count check
comfirmcases2_subjectID_check <- comfirmcases2 %>%
  group_by(subject.id) %>%
  count()

## What if we merge withdpo with comfirmcases2 directly, instead of involving the allpcr df? To do this, we will comment out the below merge that merges the allpcr and comfirmcases2 dfs.

# pcr_body_temp_symptoms <- allpcr %>%
#   inner_join(comfirmcases2, by = c('subject.id', 'sample.id', 'type.inf')) %>%
#   select(-date_on_sx)
# pcr_body_temp_symptoms$date_on_sx <- as.Date(pcr_body_temp_symptoms$date_on_sx)

finaldata <- withdpo %>%
  inner_join(comfirmcases2, by = c("sample.id", "subject.id", "type.inf", "date.visit")) %>%
  select(-date_on_sx.y) %>%
  mutate(date_on_sx = date_on_sx.x) %>%
  select(-date_on_sx.x)

# Subject ID count check on the finaldata df
finaldata_subjectID_check <- finaldata %>%
  group_by(subject.id) %>%
  count()
print(nrow(finaldata_subjectID_check))

# finaldata <- inner_join(withdpo, pcr_body_temp_symptoms, by = c("subject.id", "date.visit", "type", "type.inf", "sample.type", "final.copies", "Experiment"))

# Add body temp variable with the clinical_umd_1 df
finaldata2 <- finaldata %>%
  left_join(clinical_umd_1, by = c('date.visit', 'subject.id'))

# Add asthma variable with the clinical_umd_3 df
finaldata3 <- finaldata2 %>%
  inner_join(clinical_umd_3, by = 'subject.id')

enrolled <- readRDS("Curated Data/Cleaned Data/EMIT_subtypes_enrolled.RDS")
enrolledcase <- enrolled %>% 
  select(subject.id)
# So far this 'enrolled' object doesn't is not incoporated at all.
# It shows the list of 178 enrolled study participants

####**************G2 LOG DATA****************####
g2_in_file <- 'UMD_Raw_Data/GII/EMITGIILogUMD2013.csv'
g2_log <- read.csv(g2_in_file)

##****************** Date Entry Error Correction
print(select(filter(g2_log, subject_id == 284), subject_id, redcap_event_name, start_dt))

# Subject_id 284 g2 collection_2_arm_1 was entered as 2013-03-17 but baseline was on 2013-02-16 and collection_3 was 2013-02-18."
# Therefore recode collection_2 date to February from March (i.e. to 2013-02-17).

g2_log$start_dt[which(g2_log$subject_id == 284 & g2_log$start_dt == '2013-03-17', arr.ind = TRUE)] <- '2013-02-17'
g2_log <- g2_log %>% 
  filter(!(subject_id == 81 & redcap_event_name == 'collection_2_arm_1'))

g2_log1 <- g2_log %>% 
  select(subject_id, start_dt, g2_unit, chiller_t1, chiller_t2, chiller_t3, elbow_rh1, elbow_rh2, elbow_rh3, elbow_t1, elbow_t2, elbow_t3, cond_tin1, cond_tin2, cond_tin3, cond_tout1, cond_tout2, cond_tout3)

g2_log1$chiller.t <- (g2_log1$chiller_t1 + g2_log1$chiller_t2 + g2_log1$chiller_t3) / 3
g2_log1$elbow_rh <- (g2_log1$elbow_rh1 + g2_log1$elbow_rh2 + g2_log1$elbow_rh3) / 3
g2_log1$elbow_t <- (g2_log1$elbow_t1 + g2_log1$elbow_t2 + g2_log1$elbow_t3) / 3
g2_log1$cond_tin <- (g2_log1$cond_tin1 + g2_log1$cond_tin2 + g2_log1$cond_tin3) / 3
g2_log1$cond_tout <- (g2_log1$cond_tout1 + g2_log1$cond_tout2 + g2_log1$cond_tout3) / 3

g2_log1 <- g2_log1 %>% 
  select(subject_id, start_dt, g2_unit,chiller.t, elbow_rh,elbow_t, cond_tin,cond_tout) %>%
  rename(subject.id = subject_id, date.visit = start_dt)
g2_log1$date.visit <- as.Date(g2_log1$date.visit)

# Subject ID count check for g1_log1
g2_log1_subject_ID_check <- g2_log1 %>%
  group_by(subject.id) %>%
  count()

# Merge and clean
# Add the g2 log data variables with the g2_log1 df
finaldata4 <- finaldata3 %>%
  inner_join(g2_log1, by = c('subject.id', 'date.visit'))

# The below line of code seems to be incorrect and causes issues downstream.
# This is because this line of script only takes into consideration H3N2, Unsubtypable A, and Pandemic H1, but it doesn't take into consideration the instances where there are 2 types! (like with subject)
# finaldata4$typeAB <- ifelse(finaldata3$type.inf =='H3N2' | finaldata3$type.inf == 'Unsubtypable A' | finaldata3$type.inf == 'Pandemic H1', "A", "B")
# finaldata4$typeAB[finaldata3$subject.id == 55] <- 'A'
# finaldata4$typeAB[finaldata3$subject.id == 230] <- 'A'

# SubjectID count check for finaldata4
finaldata4_subjectID_check <- finaldata4 %>%
  group_by(subject.id) %>%
  count()

VSAS <- read.csv("UMD_Raw_Data/REDCAP/vacine_smoker_antiviral_sex.csv")

# Add sex, fluvac_cur, antiviral_24h, Smoker variables from the VSAS df 
finaldata142 <- finaldata4 %>%
  inner_join(VSAS, by = 'subject.id')

# SubjectID count check for finaldata142
finaldata142_subjectID_check <- finaldata142 %>%
  group_by(subject.id) %>%
  count()
print(nrow(finaldata142_subjectID_check))

saveRDS(finaldata142, "Curated Data/Cleaned Data/finaldata142.RDS")                         

totalsamples <- finaldata142 %>% 
  distinct(subject.id, sample.id, sample.type, .keep_all = TRUE)
totalsamples_subjectID_check <- totalsamples %>%
  group_by(subject.id) %>%
  count()

totalsubjects <- finaldata142 %>% 
  distinct(subject.id, .keep_all = TRUE)

sex <- finaldata142 %>% 
  distinct(subject.id, sex, .keep_all = TRUE) %>% 
  filter(sex == 1)

flushot <- finaldata142 %>% 
  distinct(subject.id, fluvac_cur, .keep_all = TRUE) %>% 
  filter(fluvac_cur == 1)

asthma <- finaldata142 %>% 
  distinct(subject.id, asthma, .keep_all = TRUE) %>% 
  filter(asthma == 1)

smoker <- finaldata142 %>% 
  distinct(subject.id, Smoker, .keep_all = TRUE) %>% 
  filter(Smoker == 1)

antiviral <- finaldata142 %>% 
  distinct(subject.id, anitviral_24h, .keep_all = TRUE) %>%
  filter(anitviral_24h == 1)

npsubgroup <- finaldata142 %>% 
  filter(sample.type == 'Nasopharyngeal swab') %>% 
  select(Experiment, subject.id, sample.id, type, type.inf, sample.type, final.copies, date.visit, cough_number, sneeze_number, g2.run, dpo, upper_sym, lower_sym, systemic_sym, body_temp, asthma, g2_unit, chiller.t, elbow_rh, elbow_t, cond_tin, cond_tout, sex, fluvac_cur, anitviral_24h, Smoker) %>%
  arrange(subject.id, date.visit)

npsubgroup_subjectID_check <- npsubgroup %>%
  group_by(subject.id) %>%
  count()

count1 <- npsubgroup %>% 
  distinct(subject.id, .keep_all = TRUE)

count2 <- npsubgroup %>% 
  group_by(sample.id) %>% 
  summarise(n = n())

count3 <- count2 %>% 
  filter(n > 2)

finesubgroup <- finaldata142 %>% 
  filter(sample.type == 'GII condensate NO mask')

symptom <- finesubgroup %>% 
  select(subject.id, sample.id, dpo, systemic_sym, upper_sym, lower_sym, nose_run, nose_stuf, sneeze, throat_sr, earache, chest_tight, sob, cough, malaise, headache, mj_ache, sw_fever_chill, lymph_node) %>% 
  distinct(subject.id, sample.id, dpo, systemic_sym, upper_sym, lower_sym, nose_run, nose_stuf, sneeze, throat_sr, earache, chest_tight, sob, cough, malaise, headache, mj_ache, sw_fever_chill, lymph_node)

# write.csv(symptom, "C:/Users/Jing/Desktop/symptomscoreupdate.csv") 

finesubgroup <- finesubgroup %>% 
  select(Experiment, subject.id, sample.id, type, type.inf, sample.type, final.copies, date.visit, cough_number, sneeze_number, g2.run, dpo, upper_sym, lower_sym, systemic_sym, body_temp, asthma, g2_unit, chiller.t, elbow_rh, elbow_t, cond_tin, cond_tout, sex, fluvac_cur, anitviral_24h, Smoker)

count4 <- finesubgroup %>% 
  distinct(subject.id)
count5 <- finesubgroup %>% 
  group_by(sample.id) %>% 
  summarise(n = n())
count6 <- count5 %>% 
  filter(n > 2)

coarsesubgroup <- finaldata142 %>% 
  filter(sample.type == 'Impactor 5 um NO mask') %>% 
  select(Experiment, subject.id, sample.id, type, type.inf, sample.type, final.copies, date.visit, cough_number, sneeze_number, g2.run, dpo, upper_sym, lower_sym, systemic_sym, body_temp, asthma, g2_unit, chiller.t, elbow_rh, elbow_t, cond_tin, cond_tout, sex, fluvac_cur, anitviral_24h, Smoker)

count7 <- coarsesubgroup %>% 
  distinct(subject.id)
count8 <- coarsesubgroup %>% 
  group_by(sample.id) %>% 
  summarise(n = n())
count9 <- count8 %>% 
  filter(n > 2)
# Repeated format

finesubgroup <- finesubgroup %>% 
  mutate(sampleid = gsub('^[0-9]*_', '', sample.id))
finesubgroup$NPswab <- 0
finesubgroup$Coarse <- 0
finesubgroup$Fine <- 1

coarsesubgroup <- coarsesubgroup %>% 
  mutate(sampleid = gsub('^[0-9]*_', '', sample.id))
coarsesubgroup$NPswab <- 0
coarsesubgroup$Coarse <- 1
coarsesubgroup$Fine <- 0

npsubgroup<-npsubgroup %>% 
  mutate(sampleid = gsub('^[0-9]*_', '', sample.id))
npsubgroup$NPswab <- 1
npsubgroup$Coarse <- 0
npsubgroup$Fine <- 0

finaldataset <- rbind(finesubgroup, npsubgroup, coarsesubgroup) %>% 
  ungroup 
finaldataset$final.copies <- as.numeric(finaldataset$final.copies)
finaldataset <- finaldataset %>%
  arrange(subject.id, date.visit, type, sample.type, final.copies)

# finaldataset <- finaldataset[order(finaldataset$subject.id,finaldataset$dpo), ]
# finaldataset$final.copies[is.na(finaldataset$final.copies)] <- '.'

# Need to add some manipulation to the df right here to replicate what Jing produced for final input to the SAS program

# Assess the df
finaldataset_count <- finaldataset %>%
  arrange(sample.type, subject.id, date.visit, Experiment) %>%
  group_by(sample.type, subject.id, date.visit, Experiment, type.inf) %>%
  summarise(count = n())
print(nrow(finaldataset_count))

finaldataset_subjectID_check <- finaldataset %>%
  group_by(subject.id) %>%
  count()
print(nrow(finaldataset_subjectID_check))

finaldataset_sampleID_check <- finaldataset %>%
  group_by(sample.id) %>%
  count()
print(nrow(finaldataset_sampleID_check))

# I figured out why there are only single PCR replicate results (singlicate as opposed to duplicate) - it’s because these were part of the subtyping assays and there was only enough material to run the CDC panel in singlicate (or a second extraction would have been required). For this reason, the 1st day visit NP swabs were run in singlicate (many of them). Others got put on par assays later on and have duplicates. This explains the lack of consistency in the number of pcr results reported for first visit NP swabs here. We will move forward with the dataframe that we have generated. 

# We also note that we have kept all of the data where there are multiple assays on the same sample. Tobit models will help us interpret these data, especially when there were 2 assays and one was positive while the other was negative (or there were a combination of replicates that were positive/negative)

# This df is missing the following variables (compared with Jing's output):
# centimeterheight
# kilogramweight
# BMI
# cur_asthma (how is this different than the asthma variable that is already in the df?)
# lung_sym_2pos (not sure what this is or how to compute?)
# fluvac_last2y
# bothyear
# age

# Not sure what lung_symp_2pos means, but the other variables can be taken directly from the clinical_umd df or derived from variables in the clinical_umd df. 

clinical_umd_height_weight_BMI_vax_age <- clinical_umd %>%
  filter(!is.na(age)) %>%
  mutate(height_inches = height_in + 12*(height_ft)) %>%
  mutate(height_cm = 2.54*(height_inches)) %>%
  mutate(weight_kg = 0.453592*weight) %>%
  mutate(BMI = weight_kg/((height_cm/100)^2)) %>%
  mutate(vax_bothyear = ifelse(fluvac_cur == 1 & fluvac_last2y == 1, 1, 0)) %>%
  rename(subject.id = field_subj_id, date.visit = date_visit) %>%
  select(subject.id, date.visit, height_cm, weight_kg, BMI, fluvac_last2y, vax_bothyear, fluvac_10y, age)

missing_fluvac_last2y <- clinical_umd_height_weight_BMI_vax_age %>%
  filter(is.na(fluvac_last2y))
print(nrow(missing_fluvac_last2y))

finaldataset <- finaldataset %>%
  left_join(clinical_umd_height_weight_BMI_vax_age, by = c("subject.id")) %>%
  select(-date.visit.y) %>%
  rename(date.visit = date.visit.x) %>%
  filter(subject.id != 52) %>% # It was decided by the lab that subject 52 was a false positives and should be removed
  filter(subject.id != 58) %>% # It was decided by the lab that subject 58 was a false positives and should be removed
  arrange(subject.id, date.visit, type, sample.type, Experiment, final.copies)

finaldataset_subjectID_check <- finaldataset %>%
  group_by(subject.id) %>%
  count()
print(nrow(finaldataset_subjectID_check))

write.csv(finaldataset, "Curated Data/Analytical Datasets/finaldatasetrepeatupdate.csv")

finaldataset_missing_fluvac_last2y <- finaldataset %>%
  filter(is.na(fluvac_last2y)) %>%
  distinct(subject.id)
print(nrow(finaldataset_missing_fluvac_last2y))

# Looks like 55 out of 142 subjects are missing the fluvac_last2y variable!
# There simply isn't data on this in the raw data. Is there a different raw datafile that should be used?
# Somehow in the original data that was used in the PNAS analysis, we have data for this variable on all 142 subjects.
# It looks like all of these 55 subject IDs are marked in the PNAS final dataset as having a 0 for the fluvac_last2y as opposed to an NA. I'm not sure if this is correct. 

### Creating the "all_data" df that has all the data from all of the screened and enrolled participants ####

# Need to merge together the pcr data into a definitive set called allPCRfinal
npfirst1 <- npfirst %>%
  select(subject.id, sample.id, type, Experiment, Ct, cfactor, virus.copies) %>%
  mutate(final.copies = virus.copies) %>%
  ungroup()
# The difference between copies.in and copy.num in the npfirst object is a factor of 8,000 for flu A assays and 41,100 for flu B assays.
# If we break down these factors of 8,000 (for flu A) and 41,100 (for flu B) we see that it is probably a combination of what Jing used for the conversion factor from virus particles and RNA copies (80 for flu A and 411 for flu B although it should be 250 for flu A and 272 for flu B), and the dilution factor (100 for NP swabs, except for a few where 100ul instead of 50ul was used to extract)
# So the virus.copies variable is actually the finalized pcr variable and all we needed to do in the above step was rename virus.copies as final.copies.
# However, this whole process is under review because the conversion factors between virus particles and RNA copies are not what we believe that they should be. 

total.pcr1 <- total.pcr %>%
  select(subject.id, Well.Name, type, Experiment, Ct..dRn., cfactor, virus.copies) %>%
  rename(sample.id = Well.Name)
# Unlike the first visit NP swab data discussed above in the npfirst object, the total.pcr object has undergone a different data manipulation process with respect to applying calibration factors.
# For total.pcr, the copies.in variable was multiplied by the RNA copies to virus particle conversion factors of 80 (for flu A) and 411 (for flu B) in order to get out copy.num (note - we are reviewing this because we believe that the conversion factors should really be 250 and 272 for flu A and B, respectively). Then, copy.num was multiplied by the pcr assay calibration factor to get virus.copies. Unlike in the NP swab first visit data where the dilution factor was applied as part of the process of getting from copies.in to copy.num, with the aerosols and post first visit NP swabs, the dilution factor has not yet been applied and, thus, the following code is written to get from virus.copies to a final.copies variable that has the dilution factor applied to it. 

## Need the below code (copied from above) to manipulate the pcr data from the total.pcr object in order to apply the dilution factors for the aerosol data and NP swabs that were run on visits 2 or 3 (not first visit NPS).

samples.cc_type <- samples.cc %>%
  select(sample.id, sample.type)

total.pcr1 <- total.pcr1 %>%
  left_join(enrolled, by = "subject.id") %>%
  left_join(samples.cc_type, by = "sample.id") %>%
  ungroup()

# Pick out all the PCR with A assay results
allpcrA <- total.pcr1 %>% 
  filter(type == 'A') %>%
  rename(virus.copiesA = virus.copies, typeA = type) %>%
  rename(CtA = Ct..dRn.)

# Pick out all the PCR with B assay results
allpcrB <- total.pcr1 %>% 
  filter(type == 'B') %>%
  rename(virus.copiesB = virus.copies, typeB = type) %>%
  rename(CtB = Ct..dRn.)

# Join both A and B assay, and also add sample type in the data list, assign the final RNA copies number for each sample type
allPCR <- allpcrA %>%
  full_join(allpcrB, by = c('subject.id', 'Experiment', 'type.inf', 'sample.id', 'sample.type', 'cfactor')) %>%
  arrange(subject.id) %>% 
  filter(!sample.type == 'Throat Swab')

## DATA EDITING (dilution factors different for a few samples) ##
# Flagged samples all NP swabs, dilution factor is 50
# 66_7 120_7 184_8 188_7 189_7 192_7 196_7 262_7 277_7 284_12 284_7 296_12 296_7
# Note that none of these samples were first visit NPS samples, so this data edit applied to only the data from non-first visit NPS. 

allPCR1 <- allPCR %>% 
  filter(sample.id == '66_7' | sample.id == '120_7' | sample.id == '184_8' | sample.id == '188_7' | sample.id == '189_7' | sample.id == '192_7' | 
           sample.id == '196_7' | sample.id == '262_7' | sample.id == '277_7' | sample.id == '284_12' | sample.id == '284_7' | 
           sample.id == '296_12'| sample.id=='296_7')

allPCR2 <- allPCR %>%
  anti_join(allPCR1)

allPCR3 <- allPCR1 %>% 
  mutate(final.copiesA = virus.copiesA*50, final.copiesB = virus.copiesB*50)

allPCR4 <- allPCR2 %>% 
  filter(!sample.type == 'Nasopharyngeal swab') %>% 
  mutate(final.copiesA = virus.copiesA*25, final.copiesB = virus.copiesB*25)

allPCR5 <- allPCR2 %>% 
  filter(sample.type == 'Nasopharyngeal swab') %>%
  mutate(final.copiesA = virus.copiesA*100, final.copiesB = virus.copiesB*100)

allPCRtotal <- rbind(allPCR3, allPCR4, allPCR5) %>% 
  ungroup()

allPCR.A <- allPCRtotal %>% 
  filter(typeA == 'A') %>%
  select(-CtB, -virus.copiesB, -typeB, -final.copiesB) %>%
  rename(Ct = CtA, virus.copies = virus.copiesA, type = typeA, final.copies = final.copiesA)

allPCR.B <- allPCRtotal %>% 
  filter(typeB == 'B') %>%
  select(-CtA, -virus.copiesA, -typeA, -final.copiesA) %>%
  rename(Ct = CtB, virus.copies = virus.copiesB, type = typeB, final.copies = final.copiesB)

# merge the seperated FILE FOR A and B back together, successfully make - one row for one subject.id
allPCRfinal <- rbind(allPCR.A, allPCR.B) %>% 
  ungroup() %>%
  select(-type.inf, -sample.type)

allPCRfinal_subjectID_check <- allPCRfinal %>%
  group_by(subject.id) %>%
  count()
print(nrow(allPCRfinal_subjectID_check))
# 202 subject IDs have some sort of pcr data for aerosols or post day 1 NPS

## Now, bind the allPCRfinal and npfirst1 objects together

pcr_full <- rbind(allPCRfinal, npfirst1)

pcr_full_subjectID_check <- pcr_full %>%
  group_by(subject.id) %>%
  count()
print(nrow(pcr_full_subjectID_check))
# 207 subjects

pcr_full_merge <- pcr_full %>%
  select(subject.id, sample.id, type, Experiment, final.copies)

## Now, get the samples.cc variable ready to merge.

samples.cc_date_on_sx_subjectID_check <- samples.cc %>%
  select(subject.id, sample.id, date.visit, date_on_sx) %>%
  filter(!is.na(date_on_sx)) %>%
  group_by(subject.id) %>%
  count()
print(nrow(samples.cc_date_on_sx_subjectID_check))
# We see that the samples.cc df has date_on_sx for 331 subjects (same as for the daypickfinal_full df). 

# Now let's check how many sampling instances the samples.cc has
samples.cc_sampling_instances_check <- samples.cc %>%
  distinct(subject.id, date.visit)
print(nrow(samples.cc_sampling_instances_check))
# 473 sampling instances. But what about if we exclude to only those sampling instances that are associated with g2 visits?
samples.cc_gii_sampling_instances_check <- samples.cc %>%
  filter(g2.run != 0) %>%
  distinct(subject.id, date.visit)
print(nrow(samples.cc_gii_sampling_instances_check))
# There are 276 sampling instances. 

## Now get the enrolled df ready to merge
enrolled_type.inf <- enrolled
print(nrow(enrolled_type.inf))

## Now get the daypickfinal_full df ready to merge
daypickfinal_dpo <- daypickfinal_full %>%
  select(subject.id, date.visit, dpo)
print(nrow(daypickfinal_dpo))
# Note that the date_on_sx variable in the samples.cc df is more comprehensive than the one from daypickfinal_dpo because it includes data on more subjectIDs.
# samples.cc contains 473 sampling instances, while daypickfinal_full contains 440 sampling instances.

## Now get the clinical_umd_body_temp df ready to merge
clinical_umd_body_temp <- clinical_umd_1 
print(nrow(clinical_umd_body_temp))
# Oddly enough, the clinical_umd_1 df is actually more comprehensive than the body_temp variable from the samples.cc df (this could be due to the way merges were done in the original script where samples.cc was created - this script is embedded in the current script now)

# Now get the clinical_umd_symptoms df ready to merge
clinical_umd_symptoms <- clinical_umd_2
print(nrow(clinical_umd_symptoms))

# Now get the clinical_umd_asthma df ready to merge
clinical_umd_asthma <- clinical_umd_3
print(nrow(clinical_umd_asthma))

# Now get the g2 log ready to merge
g2_log1 # Already in good form
print(nrow(g2_log1))

# Now get the vacine_smoker_antiviral_sex df ready to merge
VSAS # Already in good form
print(nrow(VSAS))

# Oddly the body temp variable from the clinical_umd object appears to be more complete than the one in the samples.cc df
# Thus we will elimnate the body_temp variable from samples.cc and use the one from clinical_umd_body_temp
samples.cc <- samples.cc %>%
  select(-body_temp)
print(nrow(samples.cc))

## Now can begin to bind all of the pieces together
all_data <- samples.cc %>%
  left_join(enrolled_type.inf, by = "subject.id") %>%
  left_join(daypickfinal_dpo, by = c("subject.id", "date.visit")) %>%
  left_join(pcr_full_merge, by = c("subject.id", "sample.id")) %>%
  left_join(clinical_umd_body_temp, by = c("subject.id", "date.visit")) %>% 
  left_join(clinical_umd_symptoms, by = c("subject.id", "date.visit")) %>%
  left_join(clinical_umd_asthma, by = "subject.id") %>%
  left_join(g2_log1, by = c("subject.id", "date.visit")) %>%
  left_join(VSAS, by = "subject.id") %>%
  arrange(subject.id, date.visit, type, sample.type, final.copies)
print(nrow(all_data))

all_data$final.copies <- as.numeric(all_data$final.copies)

# Check number of subjectIDs in all_data df
all_data_subjectID_check <- all_data %>%
  group_by(subject.id) %>%
  count()
print(nrow(all_data_subjectID_check))

all_data_gii_sample_instance_check <- all_data %>%
  filter(g2.run == 1 | g2.run == 2 | g2.run == 3) %>%
  distinct(subject.id, date.visit)
print(nrow(all_data_gii_sample_instance_check))
# Hmm, I'm getting 276 entries here but we should have 278 according to PNAS SI Table S1. 
# Looking back through the G2 log and the rest of the data, I keep seeing 276 as the correct number here.
# I'm not able to replicate the 278 number here. 

# This df is missing the following variables (compared with Jing's output):
# centimeterheight
# kilogramweight
# BMI
# cur_asthma (how is this different than the asthma variable that is already in the df?)
# lung_sym_2pos (not sure what this is or how to compute?)
# fluvac_last2y
# bothyear
# age

# Not sure what lung_symp_2pos means, but the other variables can be taken directly from the clinical_umd df or derived from variables in the clinical_umd df. 

clinical_umd_height_weight_BMI_vax_age <- clinical_umd %>%
  filter(!is.na(age)) %>%
  mutate(height_inches = height_in + 12*(height_ft)) %>%
  mutate(height_cm = 2.54*(height_inches)) %>%
  mutate(weight_kg = 0.453592*weight) %>%
  mutate(BMI = weight_kg/((height_cm/100)^2)) %>%
  mutate(vax_bothyear = ifelse(fluvac_cur == 1 & fluvac_last2y == 1, 1, 0)) %>%
  rename(subject.id = field_subj_id, date.visit = date_visit) %>%
  select(subject.id, height_cm, weight_kg, BMI, fluvac_last2y, vax_bothyear, fluvac_10y, age)

missing_fluvac_last2y <- clinical_umd_height_weight_BMI_vax_age %>%
  filter(is.na(fluvac_last2y))
print(nrow(missing_fluvac_last2y))
# Note there are 121 missing fluvac_last2y data on 121 subjects. It looks like in Jing's final dataset used in tobit regression, these NAs were turned into 0s. 

all_data <- all_data %>%
  left_join(clinical_umd_height_weight_BMI_vax_age, by = c("subject.id")) %>%
  arrange(subject.id, date.visit, type, sample.type, Experiment, final.copies)

# Note that the finaldataset only has data for NPS, fine, and coarse, so note that most of the variables for the sample.type vars levels of throat swab and anterior nasal swab are NA. 

# Write out the all_data df
write.csv(all_data, "Curated Data/Analytical Datasets/all_screened.csv")

#### Creating the all_cases dfs ####

## Now, to make the all_cases df, we take only the subjects that are in the enrolled_type.inf object and also cut to only those with g2 data
all_cases <- enrolled_type.inf %>%
  select(subject.id) %>%
  inner_join(all_data)
print(nrow(all_cases))
  
write.csv(all_cases, "Curated Data/Analytical Datasets/all_cases.csv")

all_cases_gii_samples <- enrolled_type.inf %>%
  select(subject.id) %>%
  inner_join(all_data) %>%
  filter(g2.run != 0)
print(nrow(all_cases_gii_samples))

# Check number of subjectIDs in all_cases df
all_cases_gii_samples_subjectID_check <- all_cases_gii_samples %>%
  group_by(subject.id) %>%
  count()
print(nrow(all_cases_gii_samples_subjectID_check))
# There are 178 subjects

all_cases_gii_sample_instance_check <- all_cases_gii_samples %>%
  distinct(subject.id, date.visit)
print(nrow(all_cases_gii_sample_instance_check))
# There are 276 gii sampling instances on these 178 subjects. 

# Write out the all_cases df
write.csv(all_cases_gii_samples, "Curated Data/Analytical Datasets/all_cases_gii_samples.csv")

#### Creating flu_cases dfs ####
## Now, to make the flu_cases df, we take the all the data where we have an enrolled participant (i.e., data on subtype used as a positive key for if someone enrolled and was positive for influenza virus)
flu_cases <- subtype %>%
  select(subject.id) %>%
  inner_join(all_data)
print(nrow(flu_cases))

write.csv(flu_cases, "Curated Data/Analytical Datasets/flu_cases.csv")

# We will also make another version of this that includes only the data where gii sampling instances occurred. 
flu_cases_gii_samples <- subtype %>%
  select(subject.id) %>%
  inner_join(all_data)  %>%
  filter(g2.run != 0)
print(nrow(flu_cases_gii_samples))
  
# Check number of subjectIDs in all_cases df
flu_cases_gii_samples_subjectID_check <- flu_cases_gii_samples %>%
  group_by(subject.id) %>%
  count()
print(nrow(flu_cases_gii_samples_subjectID_check))
# 158 subjectIDs

# Check the number of gii sampling instances
flu_cases_gii_samples_gii_sampling_instances_check <- flu_cases_gii_samples %>%
  distinct(subject.id, date.visit)
print(nrow(flu_cases_gii_samples_gii_sampling_instances_check))
# 250 gii sampling instances for the 150 subjects

# Let's examine this set of sampling instances that was just eliminated when we went from all_cases_gii_samples to flu_cases_gii_samples
all_cases_to_flu_cases_gii_samples <- all_cases_gii_samples %>%
  anti_join(flu_cases_gii_samples)

# Check to make sure all of these instances were negative (for the instances where there was pcr data)
all_cases_to_flu_cases_gii_samples_pcr <- all_cases_to_flu_cases_gii_samples %>%
  filter(!is.na(final.copies))
# Correct. 
# So, now let's count how many gii sampling instances are part of this all_cases_to_flu_cases_gii_samples df

all_cases_to_flu_cases_gii_samples_gii_sampling_instance_check <- all_cases_to_flu_cases_gii_samples %>%
  distinct(subject.id, date.visit)
print(nrow(all_cases_to_flu_cases_gii_samples_gii_sampling_instance_check))
# There were 26 instances here that were removed due to "not confirmed by PCR"
# Were these 26 instances just among the 20 subjects excluded or were there others gii instances excluded from some subjects outside of the grouop of 20?
# The 26 should just be from the 20 subjects excluded but let's double check
all_cases_to_flu_cases_gii_samples_gii_sampling_instance_check_subID <- all_cases_to_flu_cases_gii_samples %>%
  distinct(subject.id, date.visit) %>%
  distinct(subject.id)
print(nrow(all_cases_to_flu_cases_gii_samples_gii_sampling_instance_check_subID))
# Sure enough, the 26 gii instances came from the 20 subjects excluded.

# Write out the flu_cases_gii_samples df
write.csv(flu_cases_gii_samples, "Curated Data/Analytical Datasets/flu_cases_gii_samples.csv")

#### Getting from the 158 subjects and 250 sampling instances down to the numbers in the final PNAS dataset ####

## Exclude visits on dpo = 0
flu_cases_gii_samples_exclude_day0 <- flu_cases_gii_samples %>%
  filter(dpo != 0)

# Check subject IDs
flu_cases_gii_samples_exclude_day0_subjectID_check <- flu_cases_gii_samples_exclude_day0 %>%
  distinct(subject.id)
print(nrow(flu_cases_gii_samples_exclude_day0_subjectID_check))
# 156 subjects (2 subjects fewer than before this round of exclusion)

# Check gii sampling instances
flu_cases_gii_samples_exclude_day0_gii_sampling_instance_check <- flu_cases_gii_samples_exclude_day0 %>%
  distinct(subject.id, date.visit)
print(nrow(flu_cases_gii_samples_exclude_day0_gii_sampling_instance_check))
# 242 gii sampling instances (8 subjects fewer than before this round of exclusion)

# Who are the subjects that contributed to thes 8 gii sampling instances excluded  and how many of these instances do they account for each?
flu_cases_gii_samples_exclude_day0_gii_sampling_instance_check_subID <- flu_cases_gii_samples_gii_sampling_instances_check %>%
  anti_join(flu_cases_gii_samples_exclude_day0_gii_sampling_instance_check) %>%
  group_by(subject.id, date.visit) %>%
  count()
print(nrow(flu_cases_gii_samples_exclude_day0_gii_sampling_instance_check_subID))
# This shows that there were 8 subjects that each lost a single gii sampling instance, however 6 of these had multiple gii sampling instances and thus, not all of their data was excluded by this exclusion step.
# To check which subjects were the 2 that were completely excluded via this data exclusion step...
dpo0_subjects <- flu_cases_gii_samples_subjectID_check %>%
  anti_join(flu_cases_gii_samples_exclude_day0_subjectID_check) %>%
  count()
print(nrow(dpo0_subjects))

## Exclude visits after dpo = 3
flu_cases_gii_samples_exclude_day0_and_dpo4plus <- flu_cases_gii_samples_exclude_day0 %>%
  filter(dpo == 1 | dpo == 2 | dpo == 3)

# Check subject IDs
flu_cases_gii_samples_exclude_day0_and_dpo4plus_subjectID_check <- flu_cases_gii_samples_exclude_day0_and_dpo4plus %>%
  distinct(subject.id, .keep_all = TRUE)
print(nrow(flu_cases_gii_samples_exclude_day0_and_dpo4plus_subjectID_check))
# 149 subjects

# Check gii sampling instances
flu_cases_gii_samples_exclude_day0_and_dpo4plus_gii_sampling_instance_check <- flu_cases_gii_samples_exclude_day0_and_dpo4plus %>%
  distinct(subject.id, date.visit)
print(nrow(flu_cases_gii_samples_exclude_day0_and_dpo4plus_gii_sampling_instance_check))
# 232 gii sampling instances

# Who are the subjects that contributed to thes 10 gii sampling instances excluded  and how many of these instances do they account for each?
flu_cases_gii_samples_exclude_day0_and_dpo4plus_gii_sampling_instance_check_subID <- flu_cases_gii_samples_exclude_day0_gii_sampling_instance_check %>%
  anti_join(flu_cases_gii_samples_exclude_day0_and_dpo4plus_gii_sampling_instance_check) %>%
  group_by(subject.id, date.visit) %>%
  count()
print(nrow(flu_cases_gii_samples_exclude_day0_and_dpo4plus_gii_sampling_instance_check_subID))
# This shows that there were 10 subjects that each lost a single gii sampling instance, however 3 of these had multiple gii sampling instances and thus, not all of their data was excluded by this exclusion step.
# To check which subjects were the 7 that were completely excluded via this data exclusion step...
dpo0_and_4plus_subjects <- flu_cases_gii_samples_exclude_day0_subjectID_check %>%
  anti_join(flu_cases_gii_samples_exclude_day0_and_dpo4plus_subjectID_check) %>%
  count()
print(nrow(dpo0_and_4plus_subjects))
# To check which subjects were the other 3 that were not completely excluded via this data exclusion step...
remaining <- flu_cases_gii_samples_exclude_day0_and_dpo4plus_gii_sampling_instance_check_subID %>%
  anti_join(dpo0_and_4plus_subjects)
print(nrow(remaining))

## Exclude where cough data is missing
flu_cases_gii_samples_exclude_day0_and_dpo4plus_and_missingcough <- flu_cases_gii_samples_exclude_day0_and_dpo4plus %>%
  filter(!is.na(cough_number))

# Check subject IDs
flu_cases_gii_samples_exclude_day0_and_dpo4plus_and_missingcough_subjectID_check <- flu_cases_gii_samples_exclude_day0_and_dpo4plus_and_missingcough %>%
  distinct(subject.id)
print(nrow(flu_cases_gii_samples_exclude_day0_and_dpo4plus_and_missingcough_subjectID_check))
# 147 subjects

# Check gii sampling instances
flu_cases_gii_samples_exclude_day0_and_dpo4plus_and_missingcough_gii_sampling_instance_check <- flu_cases_gii_samples_exclude_day0_and_dpo4plus_and_missingcough %>%
  distinct(subject.id, date.visit)
print(nrow(flu_cases_gii_samples_exclude_day0_and_dpo4plus_and_missingcough_gii_sampling_instance_check))
# 225 gii sampling instances

# Who are the subjects that contributed to these 7 gii sampling instances excluded and how many of these instances do they account for each?
flu_cases_gii_samples_exclude_day0_and_dpo4plus_and_missingcough_gii_sampling_instance_check_subID <- flu_cases_gii_samples_exclude_day0_and_dpo4plus_gii_sampling_instance_check %>%
  anti_join(flu_cases_gii_samples_exclude_day0_and_dpo4plus_and_missingcough_gii_sampling_instance_check) %>%
  distinct(subject.id)
print(nrow(flu_cases_gii_samples_exclude_day0_and_dpo4plus_and_missingcough_gii_sampling_instance_check_subID))
# This shows that there were 7 subjects that each lost a single gii sampling instance, however 5 of these had multiple gii sampling instances and thus, not all of their data was excluded by this exclusion step.
# To check which subjects were the 2 that were completely excluded via this data exclusion step...
dpo0_and_4plus_cough_subjects <- flu_cases_gii_samples_exclude_day0_and_dpo4plus_subjectID_check %>%
  anti_join(flu_cases_gii_samples_exclude_day0_and_dpo4plus_and_missingcough_subjectID_check) %>%
  count()
print(nrow(dpo0_and_4plus_cough_subjects))
# To check which subjects were the other 5 that were not completely excluded via this data exclusion step...
remaining <- flu_cases_gii_samples_exclude_day0_and_dpo4plus_and_missingcough_gii_sampling_instance_check_subID %>%
  anti_join(dpo0_and_4plus_cough_subjects, by = "subject.id")
print(nrow(remaining))
## Exclude where there is incomplete PCR data

# DATA EXCLUSION DUE TO PCR & OTHER ISSUES ##
# Drop 333 only gii visit because lost coarse aerosol sample (although good data exists for the NP and the fine aerosol)
# Drop 52 because false positive (this makes 52 a negative case and thus we exclude all 2/2 gii sampling instances)
# Drop 58 because false positive (58 only had 1 gii sampling instance so this was excluded)
# Drop 182 2nd gii visit because bad interrun calibrator on the PCR (there is still a 1st gii visit for 182 so this subject is not excluded entirely)
# Drop 322 because bad interrun calibrator on the PCR (this was the only gii sampling instance so this subject is excluded entirely)
# Drop 337 because bad interrun calibrator on the PCR (this was the only gii sampling instance so this subject is excluded entirely)

flu_cases_gii_samples_exclude_day0_and_dpo4plus_and_missingcough_incompletePCR <- flu_cases_gii_samples_exclude_day0_and_dpo4plus_and_missingcough %>%
  filter(!subject.id %in% c(333, 52, 58, 322, 337)) %>%
  filter(!(subject.id == 182 & date.visit == "2013-02-08"))

## END DATA EXCLUSION DUE TO PCR & OTHER ISSUES ##

# Check subject IDs
flu_cases_gii_samples_exclude_day0_and_dpo4plus_and_missingcough_incompletePCR_subjectID_check <- flu_cases_gii_samples_exclude_day0_and_dpo4plus_and_missingcough_incompletePCR %>%
  distinct(subject.id)
print(nrow(flu_cases_gii_samples_exclude_day0_and_dpo4plus_and_missingcough_incompletePCR_subjectID_check))
# 142 subjects

# Check gii sampling instances
flu_cases_gii_samples_exclude_day0_and_dpo4plus_and_missingcough_incompletePCR_gii_sampling_instance_check <- flu_cases_gii_samples_exclude_day0_and_dpo4plus_and_missingcough_incompletePCR %>%
  filter(g2.run != 0) %>%
  distinct(subject.id, date.visit)
print(nrow(flu_cases_gii_samples_exclude_day0_and_dpo4plus_and_missingcough_incompletePCR_gii_sampling_instance_check))
# 218 gii sampling instances

# Let's make the name of this df a little better to understand and remove superfluous observations to get out what we got for the finaldataset df
PNAS_data_full <- flu_cases_gii_samples_exclude_day0_and_dpo4plus_and_missingcough_incompletePCR

## Remove the throat swab and anterior nasal swab information
## Also remove all the times that the type.inf doesn't match with the type

table(PNAS_data_full$type.inf)
table(PNAS_data_full$type)

PNAS_data_full <- PNAS_data_full %>%
  filter(sample.type != "Throat Swab") %>%
  filter(sample.type != "anterior nasal swab") %>%
  filter((type.inf == "H3N2" & type == "A") | 
           (type.inf == "Pandemic H1" & type == "A") |
           (type.inf == "Unsubtypable A" & type == "A") |
           (type.inf == "H3N2 and PH1" & type == "A") |
           (type.inf == "B" & type == "B") |
           (type.inf == "B and unsubtypable A" & (type == "B" | type == "A")) |
           (type.inf == "H3N2 and B" & (type == "B" | type == "A")))

## Compare the PNAS_data_full df with the finaldataset that was already produced.
# The only difference should be that the PNAS_data_full has some vars that aren't in finaldataset.
# Also note that the finaldataset df has indicator variables for fine, coarse, and NPS, symptoms, height, weight, BMI, vax vars, age (all variables that were added to the df at the very end of manipulation), that the PNAS_data_full doesn't have, but could be added easily. 

compare(PNAS_data_full, finaldataset)
summary(compare(PNAS_data_full, finaldataset))

# In fact, these dataset observations are identical!

# Note that the PNAS_data_full is called PNAS becasue it has additional variables, such as rapid test results, that the finaldataset (and the PNAS set that Jing used) don't have

write.csv(PNAS_data_full, "Curated Data/Analytical Datasets/PNAS_data_full.csv")

# The EMIT Main Quarantine manuscript requires fine and coarse aerosol virus RNA copy GM and SD for comparison
# To facilitate this analysis for the EMIT Main Q manuscript, we will write out a copy of this dataset to the EMIT_Data_Analysis_Jake/EMIT_Quarantine/Curated Data/Analytical Datasets directory
write.csv(PNAS_data_full, "/Users/jbueno/Box Sync/EMIT/EMIT_Data_Analysis_Jake/EMIT_Quarantine/Curated Data/Analytical Datasets/EMIT_UMD_PNAS_data_full.csv")

# The EMIT Natural Versus Artificial Inoculation manuscript requires this PNAS_data_full df as well. Thus we will write it out to the appropriate directory as well.
write.csv(PNAS_data_full, "/Users/jbueno/Box Sync/EMIT/EMIT_Data_Analysis_Jake/Natural_vs_Artificial_Infection/Analytical Datasets/EMIT_UMD_PNAS_data_full.csv")

#### Exploring the PNAS_data_full df and the flu_cases df a little more ####

## Let's do some summary analysis 

table(PNAS_data_full$sample.type)

# Number of subjects with at least one positive aerosol sample on at least one day of sampling
table(PNAS_data_full$sample.type)

PNAS_data_full_pos_aerosol <- PNAS_data_full %>%
  filter(sample.type == "GII condensate NO mask" | sample.type == "Impactor 5 um NO mask") %>%
  filter(!is.na(final.copies)) %>%
  distinct(subject.id)
print(nrow(PNAS_data_full_pos_aerosol))

# Number of subjects with at least one positive fine particle aerosol sample on at least one day of sampling
PNAS_data_full_pos_fine <- PNAS_data_full %>%
  filter(sample.type == "GII condensate NO mask") %>%
  filter(!is.na(final.copies)) %>%
  distinct(subject.id)
print(nrow(PNAS_data_full_pos_fine))

# Number of subjects with at least one positive coarse particle aerosol sample on at least one day of sampling
PNAS_data_full_pos_coarse <- PNAS_data_full %>%
  filter(sample.type == "Impactor 5 um NO mask") %>%
  filter(!is.na(final.copies)) %>%
  distinct(subject.id)
print(nrow(PNAS_data_full_pos_coarse))

# Total number of subjects
PNAS_data_full_subjects <- PNAS_data_full %>%
  distinct(subject.id)
print(nrow(PNAS_data_full_subjects))

# Fraction positive on at least one day
print(nrow(PNAS_data_full_pos_aerosol)/nrow(PNAS_data_full_subjects))

# Table with number of subjects with 1 gii instance, 2 gii instances, 3 gii instances
gii_instances <- PNAS_data_full %>%
  distinct(subject.id, g2.run) %>%
  mutate(one_instance = ifelse(g2.run == 1 & g2.run != 2 & g2.run != 3, 1, 0)) %>%
  mutate(two_instances = ifelse(g2.run == 2 & g2.run != 1 & g2.run != 3, 1, 0)) %>%
  mutate(three_instances = ifelse(g2.run == 3 & g2.run != 1 & g2.run != 2, 1, 0))
print(sum(gii_instances$one_instance))
print(sum(gii_instances$two_instances))
print(sum(gii_instances$three_instances))

table(gii_instances$g2.run)

#lapply(gii_instances, function(x) data.frame(table(x)))

# Person-visits with negative NPS
gii_person_visit_negative_NPS <- PNAS_data_full %>%
  filter(sample.type == "Nasopharyngeal swab") %>%
  filter(is.na(final.copies)) %>%
  distinct(subject.id, date.visit, .keep_all = TRUE) %>%
  select(subject.id, date.visit, type.inf)
print(gii_person_visit_negative_NPS)
# However a closer look at these 16 instances revealed that there were only sampling instances from 5 subjects where there were positive results but negative NPS data. 

# Person-visits with all negative samples
ever_positive <- PNAS_data_full %>%
  filter(!is.na(final.copies)) %>%
  distinct(subject.id, date.visit)
person_visit_negative <- PNAS_data_full %>%
  anti_join(ever_positive) %>%
  distinct(subject.id, date.visit)
print(person_visit_negative)
