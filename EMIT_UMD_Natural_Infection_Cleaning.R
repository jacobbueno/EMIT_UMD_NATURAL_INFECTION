# EMIT UMD Natural Infection Study Data Curation - cleaning raw data to produce cleaned spreadsheets
# Program Objective: Take the datasets identified as critical, clean them, and later merge to form curated one or more curated datasets
# Author: Jacob Bueno de Mesquita using material from Jing Yan and Don Milton
# Date: December 14, 2018
# Summary:

#### Load required packages and set working directory ####

library(tidyverse)
library(RcppRoll)
library(readxl)
library(knitr)
library(data.table)
library(lubridate)

setwd("/Users/jbueno/Box Sync/EMIT/EMIT_Data_Analysis_Jake")

sessionInfo() # for reproducibility

# Now pasting code from Jing Yan and Don Milton that was used in previous work on the EMIT UMD data.
# The goal here is to review their script and improve the clarity

#### **** Using Script: Jing Yan and Dr. Milton's "Merge_1-3.R-update.r" **** ####

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

#### READ in and work with the *CLINICAL DATABASE* ####

clinical_in_file <- 'EMIT_UMD_Natural_Infection/UMD_Raw_Data/REDCAP/EMITClinicalUMD2013.csv'
clinical_umd <- read.csv(clinical_in_file)

# Let's produce some summary information about this clinical_umd df

print(nrow(clinical_umd))
print(ncol(clinical_umd))

print(sum(clinical_umd$redcap_event_name=='visit_1_part_a_arm_1'))
print(sum(clinical_umd$redcap_event_name=='screen_visit_2_arm_1'))
print(sum(clinical_umd$redcap_event_name=='screen_visit_3_arm_1'))
print(sum(clinical_umd$redcap_event_name=='visit_1_part_a_arm_1') +
        sum(clinical_umd$redcap_event_name=='screen_visit_2_arm_1') +
        sum(clinical_umd$redcap_event_name=='screen_visit_3_arm_1'))

print(sum(clinical_umd$redcap_event_name=='g2_run_1_arm_1'))
print(sum(clinical_umd$redcap_event_name=='g2_run_2_arm_1'))
print(sum(clinical_umd$redcap_event_name=='g2_run_3_arm_1'))
print(sum(clinical_umd$redcap_event_name=='g2_run_1_arm_1') +
        sum(clinical_umd$redcap_event_name=='g2_run_2_arm_1') +
        sum(clinical_umd$redcap_event_name=='g2_run_3_arm_1'))

print(addmargins(with(clinical_umd, table(redcap_event_name, exclude = c()))))

# Note that one subject was enrolled twice!
print(select(filter(clinical_umd, field_subj_id == 47 | field_subj_id == 187),
             field_subj_id, date_visit, redcap_event_name, date_g2_1, rapid_flu___1, rapid_flu___2, spec_note))
# This means that there was actually one less unique invidivual than we have unique subjects ids.
# We note the above and treat subject IDs as person-illness-episodes, not persons.

# Clinical data split into screening visits and g2 runs and remerged to get one row per encounter date
clinical_min <- clinical_umd %>% 
  select(field_subj_id, redcap_event_name, date_visit, date_g2_1, rapid_flu___3, rapid_flu_loc, body_temp)

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
  select(field_subj_id, redcap_event_name, date_g2_1,g2_run) %>%
  rename(date_visit = date_g2_1)

# Find the G2 2nd and 3rd time sample collection data
clinical_g2_23 <- clinical_g2 %>%
  filter(!grepl('^g2_run_1',redcap_event_name)) %>% 
  select(field_subj_id,redcap_event_name,date_visit,g2_run)

# Merge the G2 sample collection visits (1, 2, and 3) together into the clinical_g2 df
clinical_g2 <- merge(clinical_g2_1, clinical_g2_23, 
                     c('field_subj_id','date_visit','redcap_event_name','g2_run'), all=TRUE)

sum_clinical <- merge(select(clinical_visit, -contains("date_g2_1")), clinical_g2, 
                      c('field_subj_id', 'date_visit'), all = TRUE)

sum_clinical$enrolled <- ifelse(!is.na(sum_clinical$g2_run), TRUE, FALSE)

sum_clinical <- rename(sum_clinical, visit_name = redcap_event_name.x, g2_name = redcap_event_name.y)

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
print(head(tbl_df(sum_clinical)))

#### READ in and work with the *G2 LOG DATA* ####

g2_in_file <- 'EMIT_UMD_Natural_Infection/UMD_Raw_Data/GII/EMITGIILogUMD2013.csv'
g2_log <- read.csv(g2_in_file)

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
g2_log <- filter(g2_log, !(g2_log$start_dt) == "")

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
  select(subject_id, redcap_event_name, date_visit, subj_min)

g2_log_min$g2_coll_num <- ifelse(g2_log_min$redcap_event_name == 'baseline_and_colle_arm_1', 1, 
                              ifelse(g2_log_min$redcap_event_name == 'collection_2_arm_1', 2, 3)) 

# Numbers of subjects by g2 collection event
print(ftable(addmargins(with(g2_log_min, table(redcap_event_name, g2_coll_num, exclude = c())))))
# This is an important print out of the number of G2 collection events by visit number.

g2_log_min <- g2_log_min %>% 
  select(subject_id, date_visit, g2_coll_num, subj_min) %>%
  rename(field_subj_id = subject_id)

g2_log_min$g2lm.i <- TRUE #indicator for preseence of record in g2_log_min

g2_log_min$date_visit <- as.Date(g2_log_min$date_visit)

# Check the variable names in g2_log_min
print(head(tbl_df(g2_log_min)))

#### MERGE CLINICAL AND G2-LOG DATA ####

# Here we will merge the sum_clinical df and the g2_log_min df
# These dfs were manipulated in the previous two sections of code in this script ...
# ... in preparation for this merging step.

merge1 <- merge(sum_clinical, g2_log_min, by = c('field_subj_id', 'date_visit'), all=TRUE)

# Merge 1 check the dimensions agains the source dfs
print(ftable(addmargins(with(merge1, table(clinical.i, g2lm.i, exclude=c())))))

# Select the variables of importance and their order. 
merge1 <- merge1 %>% 
  select(field_subj_id, 
         date_visit, 
         g2_coll_num, 
         enrolled, 
         visit_num, 
         g2_run, 
         clinical.i, 
         g2lm.i, 
         rapid_flu___3, 
         rapid_flu_loc, 
         body_temp)

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
  filter(indicator == 1)
merge1 <- merge1 %>% 
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
sum(merge1$visit_num == 1,na.rm = TRUE) + 
  sum(merge1$visit_num == 2,na.rm = TRUE) + 
  sum(merge1$visit_num == 3,na.rm = TRUE)

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
print(ftable(addmargins(with(merge1, table(visit_num, g2_run, g2_coll_num, exclude=c())))))

# Variables in merge1
print(head(tbl_df(merge1)))

#### READ in and work with the *FIELD SAMPLE DATABASE* ####

# Input Field Sample Data
field_db_in_file <- 'EMIT_UMD_Natural_Infection/UMD_Raw_Data/EMIT UMD Field_db/field_db.csv'
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
field.db1 <- filter(field.db1, field_subj_id != 11)
field.db1 <- filter(field.db1, field_subj_id != 28)
field.db1 <- filter(field.db1, field_subj_id != 53)
field.db1 <- filter(field.db1, field_subj_id != 73)
field.db1 <- filter(field.db1, field_subj_id != 76)

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

#### MERGE FIELD SAMPLE DATABASE WITH COMBINED CLINICAL DATABASE & G2 LOG ####

## Merge2 = merge of merge1 with field.db1 by field_subj_id and date_visit ##
merge2 <- merge(merge1, field.db1, by = c("field_subj_id", "date_visit"), all = T)

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
merge2 <- filter(merge2, !(field_subj_id == 135 & date_visit == "2013-02-05"))

# Checking merge2 rows where field.db1 not matched by rows from merge1
print(filter(merge2, is.na(merge1.i)))

# Checking: merge1 rows for subjects who had some field.db1 rows not matching merge1 rows
x <- distinct(select(filter(merge2, is.na(merge1.i)), field_subj_id))
nrow(x)

# Checking: Table of these rows
print(tbl_df(inner_join(merge1, x, by = "field_subj_id")))

# Source of data in rows of merge2 after removing 135_6.
print(ftable(addmargins(with(merge2, table(merge1.i, field.db1.i, exclude=c())))))

# Giving an indicator variable to this finalized merge2 df
merge2$merge2.i <- T

# Head of merge2
print(head(tbl_df(merge2)))

#### READ in and work with the *UMD SAMPLES DATABASE (REDCAP DATA)* ####
sample_in_file <- 'EMIT_UMD_Natural_Infection/UMD_Raw_Data/REDCAP/EMITUMDSamples2013_DATA.csv'
sample_in <- read.csv(sample_in_file, as.is = T)

sample_in$count_tech <- as.factor(sample_in$count_tech)

# Input UMD samples file (from REDCap)
sample_in_file

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

collection <- select(sample_in %>% 
                       filter(redcap_event_name == "collection_arm_1"), 
                     sample_id, field_subj_id, sample_type, date_visit )
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
passage <- select(passage, sample_id, passpos, validp)

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
print(sum(focus1$dt_count==""))
# Number of focus2 rows with no dt_count
print(sum(focus2$dt_count==""))

#### Computation of focus assay results ####

focus1$df <- 10^(ifelse(is.na(focus1$dilution_factor), 0, focus1$dilution_factor))

focus2$df <- 10^(ifelse(is.na(focus2$dilution_factor), 0, focus2$dilution_factor))

area.24 <- pi*(15.4/2)^2
area.g <- 0.64

focus1$ct_24g <- rowSums(focus1[ , c(11:20)], na.rm = T) / (10*area.g)*area.24*focus1$df/150*1000
focus1$ct_24w <- focus1$well*focus1$df/150*1000
focus1$ct_96  <- rowSums(focus1[ , c(11:13)], na.rm = T)*focus1$df/150*1000

focus1_24g <- focus1 %>%
  filter((focus1$plate_type == 1 | is.na(focus1$plate_type)) & focus1$count_meth == 1) %>% 
  select(-ct_96, -ct_24w)

focus1_24g <- focus1_24g %>%
  rename(ct = ct_24g)

focus1_24w <- focus1 %>% 
  filter((focus1$plate_type == 1 | is.na(focus1$plate_type)) & focus1$count_meth == 2) %>% 
  select(-ct_96, -ct_24g)

focus1_24w <- focus1_24w %>%
  rename(ct = ct_24w)

focus1_96  <- focus1 %>%
  filter(focus1$plate_type == 2) %>% 
  select(-ct_24w, -ct_24g)

focus1_96  <- focus1_96 %>%
  rename(ct = ct_96)

focus1_c <- arrange(rbind(focus1_96, focus1_24w, focus1_24g))

focus2$ct_24g <- rowSums(focus2[ , c(11:20)], na.rm = T) / (10*area.g)*area.24*focus2$df/150*1000

focus2$ct_24w <- focus2$well*focus2$df/150*1000

focus2$ct_96  <- rowSums(focus2[ , c(11:13)], na.rm = T)*focus2$df/150*1000

focus2_24g <- focus2 %>%
  filter((focus2$plate_type == 1 | is.na(focus2$plate_type)) & focus2$count_meth == 1) %>% 
  select(-ct_96, -ct_24w)

focus2_24g <- focus2_24g %>%
  rename(ct = ct_24g)

focus2_24w <- focus2 %>%
  filter((focus2$plate_type == 1 | is.na(focus2$plate_type)) & focus2$count_meth == 2) %>% 
  select(-ct_96, -ct_24g)

focus2_24w <- focus2_24w %>%
  rename(ct = ct_24w)

focus2_96  <- focus2 %>% 
  filter(focus2$plate_type == 2) %>% 
  select(-ct_24w, -ct_24g)

focus2_96  <- rename(focus2_96, ct = ct_96)

focus2_c <- arrange(rbind(focus2_96, focus2_24w, focus2_24g))

focus1_c <- select(focus1_c, sample_id, dt_count, count_tech, ct)

focus2_c <- select(focus2_c, sample_id, dt_count, count_tech, ct)

focus <- merge(focus1_c, focus2_c, by = "sample_id", all = T)
focus$ct <- rowMeans(cbind(focus$ct.x,focus$ct.y), na.rm = T)

missing_focus <- tbl_df(filter(focus, is.nan(ct)))
focus_allv <- focus
focus <- select(focus, sample_id, ct)

## FOCUS ASSAY RESULTS ##

# Samples listed as having a focus assay but without results
print(missing_focus)
summary(focus)

#### MERGE CULTURE RESULTS PIECE FROM FIELD SAMPLE DATABASE TO THE CUMULATIVE CLIN DB + G2 LOG + FIELD SAMPLE DB ####

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

#Number of samples where the date_visit.x (merge2) not equal date_visit.y (culture_results)
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
  filter(!(sample_id %in% c("301_3","301_5")))

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
         body_temp)

#### Write out EMIT_samples.cc.RDS file from merge of Clin DB + G2 Log + Field Sample DB ####

saveRDS(samples.cc, file = "EMIT_UMD_Natural_Infection/Curated Data/Cleaned Data/EMIT_samples.cc.RDS")

#### **** Using Script: Jing Yan and Dr. Milton's "Snippets analysis_1.r" **** ####

# Perhaps the earlier lines of script in this program address what "Snippets analysis_1.r" was mostly getting at.
# However, "Snippets analysis_1.r" provides some, perhaps useful, summary information and looks at roommates. 

# Original file info:

# By Jing Yan & Don Milton
# December 14, 2015 - December 21, 2015
# Purpose: Combine clinical redcap data with GII data from redcap and check that all subject enrolled
#          according to redcap clinical data have the appropriate number of GII records and that
#          persons not enrolled only have screening visits (up to 3) and no GII records.
#          Will also generate a list of subjects showing whether they were enrolled and how many GII
#          sessions they completed. 
#		   Will identify roommate screenings and roommates enrolled after screening as roommates. 

#### READ in and work with Clinical Database ####

clinical_umd <- read.csv(clinical_in_file)

## Check whether there was anyone with no visit 1 who has a record for having had a g2 run ##

g2_id <- unique(clinical_umd$field_subj_id[grep('^g2', clinical_umd$redcap_event_name)])

visit1_id <- unique(clinical_umd$field_subj_id[grep('^visit', clinical_umd$redcap_event_name)])

# Number of subjects with at least 1 g2 visit who also have a visit 1
print(sum(g2_id %in% visit1_id))

# Number of subjects with at least one g2 visit regardless of having a visit 1
print(length(g2_id))

# If the above two prints are the same, then everyone who gave a g2 sample had a visit 1 screening visit.

## Check to see if there are persons with visits 2 or 3 and no visit 1 ##

# Note that screening visits were only identified as screening in the redcap names for visits 2 and 3 but visit 1 was also screening

visit2_3gp <- clinical_umd %>% 
  group_by(field_subj_id, redcap_event_name) %>% 
  summarise(n = n())

visit2_3ID <- visit2_3gp$field_subj_id[grepl('^screen', visit2_3gp$redcap_event_name)]

no_visit2_3ID <- unique(visit2_3gp$field_subj_id[!grepl('^screen', visit2_3gp$redcap_event_name)])

# Subjects with Visit 2 or 3 who also have a Visit 1
print(sum(visit2_3ID %in% visit1_id))

# Number of subjects with at least one visit 2 or 3 arm regardless of having a visit 1
print(length(visit2_3ID))

## Subjects receiving second and third screening visits ##

visit2_3DT <- clinical_umd %>% 
  filter(grepl('^screen', redcap_event_name))  # all screening visit 2 and 3 records i.e. screened more than once

visit2_3DT2 <- clinical_umd %>% 
  filter(field_subj_id %in% visit2_3DT$field_subj_id) # all records for subj with a visit 2 or 3

visit2_3DT21 <- visit2_3DT2 %>% 
  filter(grepl(0, is_rmmate)) # select non roommate referrals

# nonrr_enroll_visit2_3 contains the records for people who were not roommate referrals and had visits 2 or 3.
nonrr_enroll_visit2_3 <- visit2_3DT21 %>% 
  left_join(visit2_3DT2, by = "field_subj_id") %>% 
  select(field_subj_id, redcap_event_name.y)

# Subjects 114, 127, and 250 have more than one screen visits but not referred by roomate. 
# All of them got enrolled. 114 has 2 screen visits, 127 has 3 screen visits, and 250 has 2 screen visits

# The following subjects had more than one screening visit but were not roommates.
unique(visit2_3DT21$field_subj_id)

# List of subjects and visits for repeatedly screened persons who were not roommate referrals
print(nonrr_enroll_visit2_3) 

visit2_3DT3 <- visit2_3DT2 %>% 
  filter(!grepl('^screen', redcap_event_name)) # all records except visit 2 or 3 records for subj with a visit 2 or 3

# check number of records including screening visits
# Number of subjects with a visit 2 or 3 and a visit 1 or were also renrolled
print(sum(unique(visit2_3DT2$field_subj_id) %in% unique(visit2_3DT3$field_subj_id)))

# Number of all subject with a visit 2 or 3
print(length(unique(visit2_3DT2$field_subj_id)))

## Roommates ##

# identify the roommate referal subject 
rr <- clinical_umd %>% 
  select(field_subj_id, redcap_event_name, is_rmmate, indx_id)

# rr1 gives the the list of all the subjects referred by roomate with only the visit1_arm
rr1 <- rr %>% 
  filter(rr$is_rmmate == 1)

# Total number of roommate referrals screened
print(length(rr1$field_subj_id))

# rr2 shows all the arms and all fields
rr2 <- rr1 %>% 
  left_join(rr, by = c('field_subj_id'))
rr3 <- rr2 %>% 
  select(field_subj_id, redcap_event_name.y, is_rmmate.y, indx_id.y) # limits fields from rr2

# The following roommate referred subjects were screened
print(unique(rr3$field_subj_id))

rr4 <- rr3 %>% 
  filter(grepl('^g2', redcap_event_name.y))

# The number of roommate referrals enrolled was
print(length(rr4$field_subj_id))

# The following roommate referred subjects were enrolled
print(unique(rr4$field_subj_id))
rr5 <- rr4 %>% 
  left_join(rr3, by = c('field_subj_id')) %>% 
  select(field_subj_id, redcap_event_name.y.y, is_rmmate.y.y, indx_id.y.y)

# Roommate referred subjects who were enrolled and their index ID
print(rr5)

## Method 1, find the enrolled subjects and how many times they were enrolled ##

enroll <- clinical_umd %>% 
  filter(grepl('^g2_run_1', redcap_event_name)) %>% 
  group_by(field_subj_id, redcap_event_name)

# The total number of enrolled subjects was
print(length(enroll$field_subj_id))

enroll3 <- clinical_umd %>% 
  filter(grepl('^g2_run_3', redcap_event_name)) %>% 
  group_by(field_subj_id, redcap_event_name)

# The total number of  subjects that have three GII tests was
print(length(enroll3$field_subj_id))

enroll2 <- clinical_umd %>% 
  filter(grepl('^g2_run_2', redcap_event_name)) %>% 
  group_by(field_subj_id, redcap_event_name)

# The total number of  subjects that have two GII tests was
print(length(enroll2$field_subj_id) - length(enroll3$field_subj_id))

# The total number of  subjects that have only one GII tests was
print(length(enroll$field_subj_id) - length(enroll2$field_subj_id))

# Subject 69 has 2 GII tests from redcap clinical record, but only have one GII record in the GII log on redcap, need further investigation. 

## More data exploration and manipulation ##

# remove the screen arm from the data, leave the field_subj_id and redcap_event_name and mark n = 1 for each redcap_event_name
tab0 <- clinical_umd %>% 
  filter(!grepl('^screen', redcap_event_name)) %>% 
  group_by(field_subj_id, redcap_event_name) %>% 
  summarise(n = n())

# tab1: remove screen_arm and left with visit_1 and g2_arm, and count the number of visit_1 and g2_arm
tab1 <- clinical_umd %>% 
  filter(!grepl('^screen', redcap_event_name)) %>% 
  group_by(field_subj_id, date_enroll) %>% 
  summarise(n = n())

# tab1test, remove screen_arm and only picked the field_sub_id and count the number of appearance(each subject should have at least 1 )
tab1test <- clinical_umd %>% 
  filter(!grepl('^screen', redcap_event_name)) %>% 
  group_by(field_subj_id) %>% 
  summarise(n = n())

# tab1new, add a colunm called enroll, and for n = 2 or > 2, means the subject need to have at least one time g2, so it is enrolled(1)
# otherwise it is not enrolled(0). also add a new colunm called GII_time, if if n = 2, means one g2_arm, so GII_time = n-1 = 1, apply 
# the same method for 2 and 3 GII_times
tab1new <- tab1test %>% 
  mutate(enroll = ifelse(n >= 2, 1, 0), 
         GII_time = ifelse(n >= 2, n-1, 0))

# tab1new1 sorted data for both enrolled and unenrolled subjects
tab1new1 <- tab1new %>% 
  select(field_subj_id, enroll, GII_time)
names(tab1new1)[1] <- "subject_id"

# get the enrolldate correspond with the field_subj_id
enrollDate <- clinical_umd %>% 
  filter(grepl('^visit', redcap_event_name) & date_enroll != '') %>% 
  select(field_subj_id, date_enroll)

# Number of duplicated enrolldate for the same subject id
sum(duplicated(enrollDate$field_subj_id))

# merge the enrolldate file with the ta1new which contained enroll, GII_time, by subject_id
tab1link <- tab1new %>% 
  left_join(enrollDate, by = 'field_subj_id') %>% 
  select(field_subj_id, enroll, GII_time, date_enroll)

# check if tab1link has the same number of rows as tab1new,since the row numbers are the same, and there is no missing data for date_enroll in tab1link
# so the enrolldate match with the enrollcheck and GII_times
nrow(tab1link)
nrow(tab1new) 
sum(tab1link$date_enroll == '')
names(tab1link)[1] <- "subject_id"

# subject which are enrolled
tab2link <- tab1link %>% 
  filter(tab1link$enroll == 1)

# subject came in only for screening visits
tab3link <- tab1link %>% 
  filter(tab1link$enroll == 0)
names(tab3link)[4] <- 'first_visit_date'

#### READ in and work with the G2 Log Data ####

g1 <- read.csv(g2_in_file)

# Sujbect 81 in GII file appear to have a collection_2_arm 1 but it actually is a one time GII subject
sum((g1$start_dt) == "")

# Remove the second arm of subject 81 (REVISE TO ID ROWS WITHOUT START DATE)
g2 <- filter(g1, !(g1$start_dt) == "")

# Pick the subject id and start date and order them
g3 <- g2 %>% 
  select(subject_id, start_dt) %>% 
  arrange(subject_id, start_dt)
names(g3)[2] <- "perform_date"

# Remove the duplicated and just save the first come in date
g4 <- g3[!duplicated(g2$subject_id), ]

# Group the original subject and number them
g5 <- g3 %>% 
  group_by(subject_id) %>% 
  summarise(n = n())

# g6 tells us all many times the subject enrolled for gii study
g6 <- g5 %>% 
  mutate(GII_time = ifelse(n >= 1, n, 0))

# Merge g4 (first come in date) wit the GII_time
g7 <- g4 %>% 
  left_join(g6, by = 'subject_id') %>% 
  select(subject_id, GII_time, perform_date)

# Change the start_dt as date_enroll
names(g7)[3] <- "date_enroll"

#### Producing Enrollment Summary from Clinical Database and G2 Log Data ####

# Merge the subject id with enroll time date_enroll, sample perform date, only the enrolled the subject
g8 <- g3 %>% 
  left_join(g7, by = 'subject_id') %>% 
  select(subject_id, GII_time, date_enroll, perform_date)

# Check the gii subject list with the enrolled list from culture study
m <- g7 %>% 
  left_join(tab1link, by = c('subject_id', 'GII_time', 'date_enroll'))
sum(is.na(m$GII_time))

m2 <- tab1link %>% 
  left_join(g8, by = 'subject_id', 'GII_time') %>% 
  select(subject_id, enroll, perform_date)

## Write out the enrollment summary ##
write.csv(m2, "EMIT_UMD_Natural_Infection/Curated Data/Cleaned Data/Enrollment_Summary.csv")

#### **** Using Script: Jing Yan's "field database with redcap culture.R" **** #### 

## Original information:

# "field data base check with redcap culture.R"
# by Jing Yan
# December 16, 2015
# Purpose:check if all the redcap culture sample id match with the field id   

#### READ in and work with the *FIELD SAMPLE DATABASE* ####

a <- read.csv(field_db_in_file, as.is = T)

names(a)
a1 <- a %>% 
  select(SUBJECT_IDENTIFIER, SAMPLE_ID, COLLECTION_DT, TYPE_NAME)
names(a1)[2] = "sample_id"
names(a1)[3] = "Dates_a"
names(a1)[4] = "Sample.Type"

# Add a new column as newdate which the same as collection_dt but with format m/d/y
a2 <- a1 %>% 
  mutate(newdate = as.Date(Dates_a, format = '%m/%d/%Y'))

#### READ in and work with the *UMD SAMPLES DATABASE (REDCAP DATA)* ####

# Read in redcap_culture data
b <- read.csv(sample_in_file, as.is = T)

#add a new colunm named as new date which is the same with date of sample collection
b2 <- b %>% 
  mutate(newdate = as.Date(dt_visit, format = '%m/%d/%Y'))

#b3 proves that there is no missing date of sample collection in the redcap culture file
b3 <- b2 %>% 
  filter(dt_visit != '') %>%
  rename(Sample.Type = sample_type)

m  <- b3 %>% 
  left_join(a2, by = c('sample_id', 'newdate', 'Sample.Type'))

m1 <- m %>% 
  select(sample_id, newdate, Sample.Type, Dates_a)

m11 <- m1 %>% 
  select(sample_id, newdate, Sample.Type)

sum(is.na(m1$Dates_a))

m2 <- m1 %>% 
  filter(!grepl('.', Dates_a))

# This result suggests that all the redcap culture sample IDs and sample types and enrolled data match with the field ID data ...
# ... Only 237_6 was included in the field data base but not the redcap culture data ...
# ... 237 is an enrolled subject for 1 time, it should not have a _6 sample
# From the above code we know all the culture redcap sample ID can be found in field ID, not sure if all the field ID can be found in culture redcap

m3  <- a2 %>% 
  left_join(b3, by = c('sample_id', 'newdate', 'Sample.Type'))

m4 <- m3 %>% 
  select(sample_id, newdate, Sample.Type)

sum(is.na(m4$Date.of.sample.collection))

m5 <- m4 %>% 
  filter(!grepl('.', newdate))
# m5 are the samples that were included in the field_id but not included in the redcap culture

# There is no output for this section - rather this is part of the data checking and exploratory analysis. 
# A report can be generated from the objects in this piece of the script if desired. 

#### **** Using Script: Jing and Dr. Milton's "Subtype analysis.R" script **** ####

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

#### READ in and work with the *"2016.06.17 1st visit NP swab subtyping.csv" df* ####

part1 <- read.csv('EMIT_UMD_Natural_Infection/UMD_Raw_Data/PCR Data/PCR Results/2016.06.17 1st visit NP swab subtyping.csv', as.is = T)
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

m7.s$num = NA
m7.s$num[m7.s$type.sub.H1 == T] = 1
m7.s$num[m7.s$type.sub.H3 == T] = 2
m7.s$num[m7.s$type.sub.PH1 == T] = 3
m7.s$num[m7.s$type.B == T] = 4
m7.s$num[m7.s$type.neg == T] = 5
m7.s$num[m7.s$type.H3N2.and.B == T] = 6
m7.s$num[m7.s$type.H3N2.and.PH1 == T] = 7
m7.s$num[m7.s$type.B.and.PH1 == T] = 8
m7.s$num[m7.s$type.sub.indet == T] = 9
m7.s$num[m7.s$type.badass == T] = 10

# Number of rows without a number assigned
print(nrow(filter(m7.s, is.na(num))))
# Number of rows without a number assigned - assigned to an object
check <- m7.s %>% 
  filter(is.na(num))

# Create a subset of the data (m7.s1) with subj who have only one assay or have more than one assay but the same result each time having only one row. And, with more than one row for those subj with different results for repeat assays
# Keep only rows that are different in both result 'num' and subject name
m7.s1 = m7.s[order(m7.s$subject_id), ]

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

#### READ in and work with 2016.06.17 1st visit NP swab subtyping II.csv data ####

part2 <- read.csv('EMIT_UMD_Natural_Infection/UMD_Raw_Data/PCR Data/PCR results/2016.06.17 1st visit NP swab subtyping II.csv', as.is = T)
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
  mutate(subject_id=gsub('_[1-9]*_[0-Z]*[0-Z]*[0-Z]*[0-Z]*', '', Well.Name)) %>% 
  mutate(subject_id = gsub('nf[0-Z]*', '', subject_id)) %>% 
  mutate(subject_id = gsub('dm[0-Z]*', '', subject_id)) %>%
  mutate(types = gsub('[0-9]*_[0-9]*', '', Well.Name))
n2$subject_id <- as.numeric(n2$subject_id)
n2 <- n2 %>% 
  arrange(subject_id)

finalsubtype$type.inf[finalsubtype$subject.id == 176] <- 'H3N2' 
finalsubtype$type.inf[finalsubtype$subject.id == 335] <- 'B' 
finalsubtype$type.inf[finalsubtype$subject.id == 64] <- 'Unsubtypable A' 

# Write out this finalsubtype
saveRDS(finalsubtype, file = "EMIT_UMD_Natural_Infection/Curated Data/Cleaned Data/EMIT_subtypes.RDS")

## Based on the pcr results from GII samples or 2rd/3rd np, we have modified a few subjects' subtype

updatetype <- readRDS("R_output/negative subtype sample with positive pcr.RDS")
updatetype1 <- updatetype %>% 
  select(subject.id,type) %>% distinct(subject.id, type)
updatetype2 <- updatetype1 %>% 
  filter(type == 'A')
finalsubtype$type.inf[finalsubtype$subject.id == 105] <- 'Unsubtypable A' 
finalsubtype$type.inf[finalsubtype$subject.id == 226] <- 'Unsubtypable A' 
updatetype3 <- updatetype1 %>% 
  filter(type == 'B')
# finalsubtype$type.inf[finalsubtype$subject.id == 52] <- 'B'
# finalsubtype$type.inf[finalsubtype$subject.id == 58] <- 'B'
finalsubtype$type.inf[finalsubtype$subject.id == 223] <- 'B'
finalsubtype$type.inf[finalsubtype$subject.id == 231] <- 'B'
finalsubtype$type.inf[finalsubtype$subject.id == 327] <- 'B'
finalsubtype$type.inf[finalsubtype$subject.id == 329] <- 'B'
finalsubtype$type.inf[finalsubtype$subject.id == 365] <- 'B'

enrollcheck <- readRDS("R_output/EMIT_samples.cc.RDS") %>% 
  select(subject.id,enrolled) %>% 
  distinct(subject.id,enrolled) %>% 
  filter(enrolled == TRUE)

finalenrolltype <- semi_join(finalsubytpe, enrollcheck, by = 'subject.id')
finalenrolltype <- finalenrolltype[order(finalenrolltype$subject.id), ]
finalenrolltype <- finalenrolltype %>% 
  select(subject.id, type.inf)

saveRDS(finalenrolltype, file = paste(Out.dir, "EMIT_subtypes_enrolled.RDS", sep = ""))

finalenrollepositive <- finalenrolltype %>% 
  filter(!type.inf == 'Negative')

saveRDS(finalenrollepositive, file = paste(Out.dir, "EMIT_subtypes_enrolled_positive.RDS", sep = ""))
negative <- finalenrolltype %>% 
  filter(type.inf == 'Negative')
h3n2 <- finalenrolltype %>% filter(type.inf == 'H3N2')
B <- finalenrolltype %>% 
  filter(type.inf == 'B')
Pandemic.H1 <- finalenrolltype %>% 
  filter(type.inf == 'Pandemic H1')

