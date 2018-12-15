# EMIT UMD Natural Infection Study Data Curation - cleaning raw data to produce cleaned spreadsheets
# Program Objective: Take the datasets identified as critical, clean them, and later merge to form curated one or more curated datasets
# Author: Jacob Bueno de Mesquita
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

#### "Merge_1-3.R-update" ####

# By Jing Yan & Don Milton
# January 20-25, 2016
# Purpose: follow the data analysis plan in folder EMIT_Data_Analysis.

##
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
# Read in the clinical database

clinical_in_file <- 'EMIT_UMD_Natural_Infection/UMD_Raw_Data/REDCAP/EMITClinicalUMD2013.csv'
clinical_umd<- read.csv(clinical_in_file)

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

print(addmargins(with(clinical_umd,table(redcap_event_name,exclude = c()))))

# Note that one subject was enrolled twice!
print(select(filter(clinical_umd,field_subj_id==47|field_subj_id==187),
             field_subj_id,date_visit,redcap_event_name, date_g2_1,rapid_flu___1,rapid_flu___2,spec_note))
# This means that there was actually one less unique invidivual than we have unique subjects ids.
# Subject IDs are for person-illness-episodes, not persons.

# Clinical data split into screening visits and g2 runs and remerged to get one row per encounter date
clinical_min <- clinical_umd%>%select(field_subj_id,redcap_event_name,date_visit,date_g2_1,rapid_flu___3,rapid_flu_loc,body_temp)
clinical_visit= clinical_min%>%filter(grepl('visit',redcap_event_name))
clinical_g2 = clinical_min%>%filter(grepl('^g2_run',redcap_event_name))
clinical_visit$visit_num <- ifelse(clinical_visit$redcap_event_name=='visit_1_part_a_arm_1', 1, 
                                   ifelse(clinical_visit$redcap_event_name=='screen_visit_2_arm_1',2,3)) 
clinical_g2$g2_run <- ifelse(clinical_g2$redcap_event_name=='g2_run_1_arm_1', 1, 
                             ifelse(clinical_g2$redcap_event_name=='g2_run_2_arm_1',2,3))
clinical_g2_1 <- clinical_g2%>%filter(grepl('^g2_run_1',redcap_event_name))
clinical_g2_23 <- filter(clinical_g2,!grepl('^g2_run_1',redcap_event_name))
clinical_g2_23 <- clinical_g2_23%>%select(field_subj_id,redcap_event_name,date_visit,g2_run)

clinical_g2_1 <- select(clinical_g2_1,field_subj_id,redcap_event_name,date_g2_1,g2_run)
clinical_g2_1 <- rename(clinical_g2_1,date_visit=date_g2_1)
clinical_g2 <- merge(clinical_g2_1, clinical_g2_23,c('field_subj_id','date_visit','redcap_event_name','g2_run'),all=TRUE)
sum_clinical=merge(select(clinical_visit,-contains("date_g2_1")), clinical_g2,c('field_subj_id','date_visit'),all=TRUE)
sum_clinical$enrolled<- ifelse(!is.na(sum_clinical$g2_run),TRUE,FALSE)
sum_clinical <- rename(sum_clinical, visit_name = redcap_event_name.x, g2_name = redcap_event_name.y)
sum_clinical$g2_run <- with(sum_clinical,ifelse(is.na(g2_run),0,g2_run)) 
#sum_clinical$g2_name <- with(sum_clinical,ifelse(g2_run==0,"g2_not_run",g2_name)) 
sum_clinical$visit_num <- with(sum_clinical,ifelse(is.na(visit_num),999,visit_num))
#sum_clinical$visit_name <- with(sum_clinical,ifelse(visit_num==999,"g2_run_f/u",visit_name))
sum_clinical$clinical.i <- TRUE #indicator for presence of record in sum_clinical

cat("\nNumber of rows in summary data (total number of unique encouters): ", nrow(sum_clinical), "\n")
cat("\nSummary clinical data merged the first g2 run with the screening visit from which the participant was enrolled.\n")
cat("Therefore, the number rows is reduced by the number of initial g2 runs in the clinical data.\n")

cat("\n Tabulations for data checks \n")
print(addmargins(with(sum_clinical,table(visit_name,visit_num,exclude = c()))))
cat("\n")
print(addmargins(with(sum_clinical,table(g2_name,g2_run,exclude = c()))))
cat("\n")
print(addmargins(with(sum_clinical,table(enrolled,clinical.i,exclude = c()))))
cat("\n")
print(addmargins(with(sum_clinical,table(visit_num,g2_run,exclude = c()))))
cat("\n")
cat("\n", "Total number of g2 runs (sum of runs 1, 2, and 3) according to initial clinical data: ", sum(sum_clinical$g2_run>0), "\n")
cat("\n Variables in sum_clinical: \n")
sum_clinical$date_visit <- as.Date(sum_clinical$date_visit)
print(head(tbl_df(sum_clinical)))

cat("\n***#2 a,b,c,d***   \n")

cat("\n**************G2 LOG DATA****************\n")
g2_in_file <- 'InputFiles_UMD/EMITGIILogUMD2013.csv'
g2_log <- read.csv(g2_in_file)

cat("\n Input G-II log file: ", g2_in_file)
cat("\n", "Number of rows in g2_log data: ", nrow(g2_log))
cat("\n", "Number of cols in g2_log data: ",ncol(g2_log), "\n")

cat("\n****************** Date Entry Error Correction\n")
print(select(filter(g2_log,subject_id==284),subject_id,redcap_event_name, start_dt))
cat("\nSubject_id 284 g2 collection_2_arm_1 was entered as 2013-03-17 but baseline was on 2013-02-16 and collection_3 was 2013-02-18.",
    "Therefore recode collection_2 date to February from March (i.e. to 2013-02-17).")
g2_log$start_dt[which(g2_log$subject_id==284 & g2_log$start_dt=='2013-03-17',arr.ind=TRUE)] <- '2013-02-17'
cat("\n*****************\n")


cat("\n", "Number of subject(s) without a start_dt: ", sum((g2_log$start_dt)==""))
id <- as.integer(select(filter(g2_log,start_dt==""),subject_id)) 
cat("\n", "G-II log for cases without start date:  ")
print(tbl_df(select(filter(g2_log,subject_id==id),subject_id,redcap_event_name,start_dt,g2_unit,operator,chiller_t1,subj_min)))
cat("\n", "Clinical data for cases without start date:  ")
print(tbl_df(filter(sum_clinical,field_subj_id==id)))

cat("\n****************** Data Entry Error Correction \nRemove empty record identified as having no start_dt, subject_id: ", 
    g2_log$subject_id[which(g2_log$start_dt=="", arr.ind=TRUE)])
g2_log<- filter(g2_log,!(g2_log$start_dt)=="")
cat("\n*****************\n")

cat("\n", "Number of rows in g2_log data: ", nrow(g2_log))
cat("\n", "Number of cols in g2_log data: ",ncol(g2_log), "\n")

#  cat("\n", "Number of subjects with a coll_arm_1:",sum(g2_log$redcap_event_name=='baseline_and_colle_arm_1'))
#  cat("\n", "Number of subjects with a coll_2 arm:",sum(g2_log$redcap_event_name=='collection_2_arm_1'))
#  cat("\n", "Number of subjects with a coll_3 arm:",sum(g2_log$redcap_event_name=='collection_3_arm_1'),"\n")
g2_log<- rename(g2_log,date_visit=start_dt)
g2_log_min=g2_log%>%select(subject_id,redcap_event_name,date_visit,subj_min)
g2_log_min$g2_coll_num=ifelse(g2_log_min$redcap_event_name=='baseline_and_colle_arm_1', 1, 
                              ifelse(g2_log_min$redcap_event_name=='collection_2_arm_1',2,3)) 

cat("\n Numbers of subjects by g2 collection event \n")
print(ftable(addmargins(with(g2_log_min,table(redcap_event_name,g2_coll_num,exclude = c())))))

g2_log_min=g2_log_min%>%select(subject_id,date_visit,g2_coll_num,subj_min)
g2_log_min<- rename(g2_log_min,field_subj_id=subject_id)
g2_log_min$g2lm.i <- TRUE #indicator for preseence of record in g2_log_min

g2_log_min$date_visit<- as.Date(g2_log_min$date_visit)

cat("\n Variables in g2_log_min: \n")
print(head(tbl_df(g2_log_min)))


cat("\n\n***#3 a,b,c,d*** \n")

cat("\n****************MERGE CLINICAL AND G2-LOG DATA********************  \n")
merge1=merge(sum_clinical,g2_log_min,c('field_subj_id','date_visit'),all=TRUE)
cat("\n Merge 1 check for record sources \n")
print(ftable(addmargins(with(merge1,table(clinical.i,g2lm.i,exclude=c())))))

merge1=merge1%>%select(field_subj_id,date_visit,g2_coll_num,enrolled,visit_num,g2_run,clinical.i,g2lm.i,rapid_flu___3,rapid_flu_loc,body_temp)

cat("\n", "Do the g2_coll_num match with the g2_run? They should match: ", identical(merge1$g2_coll_num, merge1$g2_run)) 

# If there is no g2 run then the collection number and run number are set to zero here. 
merge1$g2_coll_num[is.na(merge1$g2_coll_num)]<-0
merge1$g2_run[is.na(merge1$g2_run)]<-0
merge1$indicator<-ifelse(merge1$g2_coll_num==merge1$g2_run, 1,0) 
i <- which(merge1$indicator==0,arr.ind=TRUE)
cat("\n", "Since they don't match, find the unmatched subject, and shows which row and column the subject(s) located: ",i)
cat("\n", "Subject that does not match field_subj_id =", merge1$field_subj_id[i], "\n")

cat("\n****************** Data Editing")
cat("\n", "Subject 69 (field_subj_id =", merge1$field_subj_id[which(merge1$indicator==0,arr.ind=TRUE)], ") had a second g2 visit in clinical data ",
    "but did not provide sample \n -- see REDCap comments for details -- removed from final analysis data sets.\n")
merge1 <- merge1%>%filter(indicator==1)
merge1=merge1%>%select(-indicator)
cat("*********************************\n")

cat("\n Recheck after delection: Merge 1 record sources \n")
print(ftable(addmargins(with(merge1,table(clinical.i,g2lm.i,exclude=c())))))

cat("\n", "Number of rows in edited merged data: ", nrow(merge1), "\n")
cat("\n", "Number of subjects with a first visit:",sum(merge1$visit_num==1,na.rm = TRUE))
cat("\n", "Number of subjects with a second screen:",sum(merge1$visit_num==2,na.rm = TRUE))
cat("\n", "Number of subjects with a third screen:",sum(merge1$visit_num==3,na.rm = TRUE))
cat("\n", "Number of unique subjects:",length(unique(merge1$field_subj_id)))
cat("\n", "Number of screening visits: ", sum(merge1$visit_num==1,na.rm = TRUE)+sum(merge1$visit_num==2,na.rm = TRUE)
    +sum(merge1$visit_num==3,na.rm = TRUE), "\n")

cat("\n", "Number of subjects with a 1st g2 run; g2_num: ",sum(merge1$g2_run==1), " g2_coll_num: ", sum(merge1$g2_coll_num==1))
cat("\n", "Number of subjects with a 2nd g2 run; g2_num: ",sum(merge1$g2_run==2), " g2_coll_num: ", sum(merge1$g2_coll_num==2))
cat("\n", "Number of subjects with a 3rd g2 run; g2_num: ",sum(merge1$g2_run==3), " g2_coll_num: ", sum(merge1$g2_coll_num==3))
cat("\n", "Total number of g2 runs; based on g2_num: ", sum(!merge1$g2_run==0), 
    " based on g2_coll_num: ", sum(!merge1$g2_coll_num==0), "\n")
cat("\n", "Total number of screenings without a g2 run whether or not later enrolled. ")
cat("  Based on g2_num: ", sum(merge1$g2_run==0), " based on g2_coll_num: ", sum(merge1$g2_coll_num==0))
t1 <- sum(!merge1$g2_run==0)+sum(merge1$g2_run==0)
t2 <- sum(!merge1$g2_coll_num==0)+sum(merge1$g2_coll_num==0)
cat("\n", "Total number of encounters based on g2_run: ", t1, "based on g2_coll_num: ", t2, "\n")
merge1$merge1.i <- T #indicator for record in merge1

cat("\n Cross tab tables of data in merge1 \n")
print(ftable(addmargins(with(merge1,table(visit_num,g2_run,g2_coll_num,exclude=c())))))
cat("\n Variables in merge1: \n")
print(head(tbl_df(merge1)))

cat("\n\n***#4. a,b,c***\n")

cat("\n", "*******************FIELD SAMPLE DATABASE******************** \n")
field_db_in_file <- 'InputFiles_UMD/field_db.csv'
field.db<- read.csv(field_db_in_file,as.is=T)
cat("\n Input Field Sample Data: \n")
print(head(tbl_df(field.db)))
field.db1 <- field.db%>%select(SUBJECT_IDENTIFIER, SAMPLE_ID,COLLECTION_DT ,TYPE_NAME )
field.db1<- rename(field.db1,field_subj_id=SUBJECT_IDENTIFIER)
field.db1<- rename(field.db1,sample_id=SAMPLE_ID)
field.db1<- rename(field.db1,date_visit=COLLECTION_DT)
field.db1<- rename(field.db1,sample_type=TYPE_NAME)
field.db1$field.db1.i <-TRUE #indicator that data is in field.db1

cat("\n Input field sample field.db file: ", field_db_in_file)
cat("\n", "number of rows in field sample database:",nrow(field.db1))
cat("\n", "number of columns in database (selected columns):",ncol(field.db1),"\n")

cat("\n Tabluation of number of rows by sample type in field database: \n")
print(addmargins(with(field.db1,table(sample_type,field.db1.i,exclude = c()))))

cat("\n", "Number of rows that have a missing sample_type:", sum(field.db1$sample_type==""))

# commented out because some of these samples turned out to exist in other data streams 
#   cat("\n****************** Data Editing")
#   cat("\n", "Rows with missing sample type deleted.")
#   field.db1<- filter(field.db1,!(field.db1$sample_type)=="")
#   cat("\n*****************\n")

#  print(addmargins(with(field.db1,table(sample_type,field.db1.i,exclude = c()))))

cat("\n", "Number of rows that have a missing date_visit:",  sum(field.db1$date_visit==""))
cat("\n", "Number of rows that have a NA for date_visit:",  sum(is.na(field.db1$date_visit)))
cat("\n", "Number of subjects that have a missing sample_id:", sum(field.db1$sample_id==""))
cat("\n", "Number of subjects that have a NA for sample_id:", sum(is.na(field.db1$sample_id)),"\n")

field.db1$date_visit <- as.Date(field.db1$date_visit,format="%m/%d/%Y")


cat("\n****************** Data Editing\n")
cat("\n Field_subj_id 225 was moved to field_subj_id 250 in clinical database because second screening visit was ")
cat("\n erroneously given a new ID number. However the samples are still shown in the field database as 225. ")
cat("\n Therefore, I am recoding the subject id to 250 but leaving the sample_id as 225_x -- at least for now.")

field.db1$field_subj_id <- with(field.db1, ifelse(field_subj_id==225,250,field_subj_id))
print(tbl_df(filter(field.db1,field_subj_id==250)))

cat("\n ***Recode date_visit subj 10 sample_id=10_6\n")
cat("\n   orginal data:\n")
print(tbl_df(filter(field.db1,field_subj_id==10)))
field.db1$date_visit <- with(field.db1, ifelse(field_subj_id==10 & date_visit==as.Date("2012-12-07"),as.Date("2012-12-05"),date_visit))
field.db1$date_visit <- as.Date(field.db1$date_visit,format = "%Y-%m-%d", origin = "1970-01-01")
cat("\n   recoded data:\n")
print(tbl_df(filter(field.db1,field_subj_id==10)))

cat("\n ***Recode date_visit subj 12 sample_id=12_6\n")
cat("\n   orginal data:\n")
print(tbl_df(filter(field.db1,field_subj_id==12)))
field.db1$date_visit <- with(field.db1, ifelse(field_subj_id==12 & date_visit==as.Date("2013-02-08"),as.Date("2012-12-10"),date_visit))
field.db1$date_visit <- as.Date(field.db1$date_visit,format = "%Y-%m-%d", origin = "1970-01-01")
cat("\n   recoded data:\n")
print(tbl_df(filter(field.db1,field_subj_id==12)))

cat("\n ***There was no subject 11, 28, 53, 73, or 76samples generated in error.")
cat("\n Delete samples for non-existant subjects.\n")
field.db1 <- filter(field.db1,field_subj_id!=11)
field.db1 <- filter(field.db1,field_subj_id!=28)
field.db1 <- filter(field.db1,field_subj_id!=53)
field.db1 <- filter(field.db1,field_subj_id!=73)
field.db1 <- filter(field.db1,field_subj_id!=76)

cat("\n ***There was no second visit for subj 30, samples generated in error.")
cat("\n Delete samples for subject 30 on 2012-12-18. Original data: \n")
print(tbl_df(filter(field.db1,field_subj_id==30)))
field.db1 <- filter(field.db1,!(field_subj_id==30 & date_visit==as.Date("2012-12-18"))) 
cat("\n   corrected data subject 30:\n")
print(tbl_df(filter(field.db1,field_subj_id==30)))

cat("\n ***There was no second visit for subj 120 on 2013-02-08 and sample 120_12 is not in REDCap sample database.")
print(tbl_df(filter(field.db1,field_subj_id==120)))
cat("\n    Sample 120_12 deleted.")
field.db1 <- filter(field.db1,sample_id!="120_12")

cat("\n ***Samples (NP only) from erroneous second g2 visit for subject 69 (see above) deleted.\n")
field.db1 <- filter(field.db1,sample_id!="69_6" & sample_id!="69_7")

cat("\nEND****************** Data Editing Field Sample Database \n")

cat("\n", "Number of columns in EDITED filed sample database (selected columns):",ncol(field.db1),"\n")
cat("\n Variable names in field.db1: \n")
print(head(tbl_df(field.db1)))


cat("\n***#5 a-b Merge2***\n")
cat("\n ********************* Merge2 = merge of merge1 with field.db1 by field_subj_id and date_visit ******\n")
merge2 <- merge(merge1, field.db1, by=c("field_subj_id", "date_visit"), all=T)
cat("\n Head of merge2 \n")
print(head(tbl_df(merge2)))
cat("\n Source of data in rows of merge2 \n")
print(ftable(addmargins(with(merge2,table(merge1.i,field.db1.i,exclude=c())))))

cat("Remove all samples that were assigned in the field db but not collected from the unenrolled subjects.")
merge2<- filter(merge2, !(enrolled==F & sample_type %in% c("GII condensate NO mask","Throat Swab","Impactor 5 um NO mask","anterior nasal swab")))

cat("\n Source of data in rows of merge2 after removing extraneous samples.\n")
print(ftable(addmargins(with(merge2,table(merge1.i,field.db1.i,exclude=c())))))

cat("\n Merge2: rows where merge1 was not matched by rows from field.db1\n")
print(filter(merge2,is.na(field.db1.i)))

cat("\n All rows in merge2 for subjects that have some merge1 rows not matching field.db1\n")
x <- distinct(select(filter(merge2,is.na(field.db1.i)),field_subj_id))
print(inner_join(merge2,x,by="field_subj_id"))

cat("\n All rows in field.db1 for subjects that had some merge1 rows not matching field.db1\n")
print(inner_join(field.db1,x,by="field_subj_id"))

cat("\n Subject 135 returned for a second screening visit but was never enrolled. There is nothing ")
cat("\n in the REDCap clinical record to explain why no samples in the field DB were associated with the second visit.")
cat("\n There are also no samples in the REDCap sample database for subject 135 except 135_1.")
cat("\n field_subj_id==135 & date_visit==2013-02-05 deleted.\n")
merge2 <- filter(merge2, !(field_subj_id==135 & date_visit=="2013-02-05"))

cat("\n Merge2 rows where field.db1 not matched by rows from merge1\n")
print(filter(merge2,is.na(merge1.i)))

cat("\n Merge1 rows for subjects who had some field.db1 rows not matching merge1 rows: nrow= ")
x <- distinct(select(filter(merge2,is.na(merge1.i)),field_subj_id))
cat(nrow(x), " Table of these rows: \n")
print(tbl_df(inner_join(merge1,x,by="field_subj_id")))

cat("\n Source of data in rows of merge2 after removing 135_6.\n")
print(ftable(addmargins(with(merge2,table(merge1.i,field.db1.i,exclude=c())))))

merge2$merge2.i <- T
cat("\n Head of merge2 \n")
print(head(tbl_df(merge2)))


cat("\n\n***# 6, a-j***\n") 

cat("\n", "*********************UMD SAMPLES DATABASE (REDCap)********************** \n")
sample_in_file <- 'InputFiles_UMD/EMITUMDSamples2013_DATA.csv'
sample_in <- read.csv(sample_in_file,as.is=T)
sample_in$count_tech <- as.factor(sample_in$count_tech)
cat("\n Input UMD samples file (from REDCap): ", sample_in_file)
cat("\n", "Number of rows: ", nrow(sample_in))
cat("\n", "Number of cols: ", ncol(sample_in),"\n")
sample_in <- rename(sample_in, date_visit=dt_visit)
sample_in$date_visit <- as.Date(sample_in$date_visit,format="%m/%d/%Y")

cat("\n Head of sample_in dataframe after renaming date_visit and setting count_tech as factor\n")
print(tbl_df(sample_in))

cat("\n", "Number of rows that have a missing sample_id:", sum(sample_in$sample_id==""))
cat("\n", "Number of rows that have a NA for sample_id:", sum(is.na(sample_in$sample_id)))

cat("\n****************** Data Editing\n")
cat("\n Samples 69_6 and 69_7 collected for a g2_run = 2 that was not performed. See above. Deleted here.")
cat("\n Also, samples 20-1 & 97-9 are typos and duplications. They are also deleted here.")
sample_in <- filter(sample_in,!(sample_id %in% c("69_6","69_7","20-1","97-9")))
sample_in$dt_stained <- with(sample_in,ifelse(sample_id=="20_1"&dt_count=="12/15/2012","12/15/2012",dt_stained))

collection <- select(sample_in%>%filter(redcap_event_name=="collection_arm_1"), sample_id, field_subj_id, sample_type, date_visit )
assay1 <- sample_in%>%filter(grepl('^assay1',redcap_event_name))
assay2 <- sample_in%>%filter(grepl('^assay2',redcap_event_name))

cat("\n", "Number of collection records: ", nrow(collection))
cat("\n", "Number of assay 1 records: ", nrow(assay1))
cat("\n", "Number of assay 2 records: ", nrow(assay2), "\n")
cat("\n", "Event name and sample_type read in from REDCap sample database: \n")
print(ftable(addmargins(with(sample_in,table(redcap_event_name,sample_type,exclude = c())))))

cat("\n", "Head of sample_in samples with no sample_type: \n")
print(head(filter(sample_in,sample_type=="")))


passage <- select(filter(sample_in, !is.na(passage_1)), sample_id, passage_id, passage_id_problem, dt_pass, pass_tech, 
                  passage_1, passage_2, dt_pass_2, passage_complete)
cat("\ ***PROBLEM OBS NEEDS EDITING ***
    Sample_id does not match Passage_id. Refer for lab review. 
    Meanwhile, will use sample_id as it is not duplicated.")
print(tbl_df(filter(passage,sample_id != passage_id)))

passage$passpos <- (passage$passage_1==2 | passage$passage_2==2) #Passage is + if either passage is +
passage$validp <- !is.na(passage$passpos)
passage<-select(passage,sample_id,passpos,validp)

cat("\n Initial look at passage assays")
cat("\n", "Number of passage assays: ", nrow(passage))
cat("\n", "Number of passage assays with missing date of passage: ", sum(passage$dt_pass==""))
cat("\n", "Number of valid passage assays: ", sum(passage$validp))
cat("\n", "Number of invalid passage assays: ", sum(!passage$validp))
cat("\n", "Number of positive passage assays: ", sum(passage$passpos,na.rm=T))
cat("\n", "Number of negative passage assays: ", sum(!passage$passpos,na.rm=T), "\n")

focus1 <- filter(assay1, !dt_count=="")[,c(1,17:37)]
focus2 <- filter(assay2, !dt_count=="")[,c(1,17:37)]

cat("\n", "Focus Assays")
cat("\n", "Samples with miss match sample_id and focus_id focus1: ", sum(!(focus1$sample_id==focus1$focus_id)))
cat("\n", "Samples with miss match sample_id and focus_id focus2: ", sum(!(focus2$sample_id==focus2$focus_id)))
cat("\n", "Number of focus1 counts: ", nrow(focus1))
cat("\n", "Number of focus2 counts: ", nrow(focus2))
cat("\n", "Number of focus1 rows with no dt_count: ", sum(focus1$dt_count==""))
cat("\n", "Number of focus2 rows with no dt_count: ", sum(focus2$dt_count==""),"\n")



#______Computation of focus assay results______________
focus1$df <- 10^(ifelse(is.na(focus1$dilution_factor),0,focus1$dilution_factor))
focus2$df <- 10^(ifelse(is.na(focus2$dilution_factor),0,focus2$dilution_factor))
area.24 <- pi*(15.4/2)^2
area.g <- 0.64


focus1$ct_24g <- rowSums(focus1[,c(11:20)],na.rm=T) / (10*area.g)*area.24*focus1$df/150*1000
focus1$ct_24w <- focus1$well*focus1$df/150*1000
focus1$ct_96  <- rowSums(focus1[,c(11:13)], na.rm=T)*focus1$df/150*1000

focus1_24g <- filter(focus1, (focus1$plate_type==1 | is.na(focus1$plate_type)) & focus1$count_meth==1) %>% select(-ct_96, -ct_24w)
focus1_24g <- rename(focus1_24g, ct=ct_24g)

focus1_24w <- filter(focus1, (focus1$plate_type==1 | is.na(focus1$plate_type)) & focus1$count_meth==2) %>% select(-ct_96, -ct_24g)
focus1_24w <- rename(focus1_24w, ct=ct_24w)

focus1_96  <- filter(focus1, focus1$plate_type==2) %>% select(-ct_24w, -ct_24g)
focus1_96  <- rename(focus1_96, ct=ct_96)

focus1_c <- arrange(rbind(focus1_96,focus1_24w,focus1_24g))


focus2$ct_24g <- rowSums(focus2[,c(11:20)],na.rm=T) / (10*area.g)*area.24*focus2$df/150*1000
focus2$ct_24w <- focus2$well*focus2$df/150*1000
focus2$ct_96  <- rowSums(focus2[,c(11:13)], na.rm=T)*focus2$df/150*1000

focus2_24g <- filter(focus2, (focus2$plate_type==1 | is.na(focus2$plate_type)) & focus2$count_meth==1) %>% select(-ct_96, -ct_24w)
focus2_24g <- rename(focus2_24g, ct=ct_24g)

focus2_24w <- filter(focus2, (focus2$plate_type==1 | is.na(focus2$plate_type)) & focus2$count_meth==2) %>% select(-ct_96, -ct_24g)
focus2_24w <- rename(focus2_24w, ct=ct_24w)

focus2_96  <- filter(focus2, focus2$plate_type==2) %>% select(-ct_24w, -ct_24g)
focus2_96  <- rename(focus2_96, ct=ct_96)

focus2_c <- arrange(rbind(focus2_96,focus2_24w,focus2_24g))


focus1_c <- select(focus1_c,sample_id, dt_count, count_tech, ct)
focus2_c <- select(focus2_c,sample_id, dt_count, count_tech, ct)
focus <- merge(focus1_c,focus2_c, by="sample_id",all=T)
focus$ct <- rowMeans(cbind(focus$ct.x,focus$ct.y),na.rm=T)
missing_focus <- tbl_df(filter(focus,is.nan(ct)))
focus_allv <- focus
focus <- select(focus, sample_id, ct)

cat("\n***FOCUS ASSAY RESULTS***\n")
cat("\n", "Samples listed as having a focus assay but without results.  ")
print(missing_focus)
summary(focus)

cat("\n***MERGE CULTURE RESULTS***\n")
culture_results <- merge(collection, passage, by="sample_id", all=T)
culture_results <- merge(culture_results, focus, by="sample_id", all=T)
cat("\n", "N rows collection: ", nrow(collection))
cat("\n", "N rows passage: ", nrow(passage))
cat("\n", "N rows focus: ", nrow(focus))
cat("\n", "N rows culture_results: ", nrow(culture_results), "\n")

culture_results$np <- ifelse(culture_results$sample_type=='Nasopharyngeal swab',T,F)
culture_results$impactor <- ifelse(culture_results$sample_type=='Impactor 5 um NO mask',T,F)
culture_results$condensate <- ifelse(culture_results$sample_type=='GII condensate NO mask',T,F)
culture_results$antnasal <- ifelse(culture_results$sample_type=='anterior nasal swab',T,F)  
culture_results$throat <- ifelse(culture_results$sample_type=='Throat Swab',T,F)  
culture_results$focus.i <- ifelse(is.na(culture_results$ct),F,T)
culture_results$passage.i <- ifelse(is.na(culture_results$validp),F,T)
culture_results$cr.i <- TRUE #Indicator for record present in culture_results

cat("\n Samples with (1) and without (0) passage assays and focus assays \n")
print(ftable(addmargins(with(culture_results,table(sample_type,passage.i,focus.i,exclude=c()))))) 
cat("\n  *Note: All samples should have a sample type; anterior nasal swabs & impactors should not have culture assays\n")

cat("\n ********************************************************************************")
cat("\n ******PROBLEMATIC SAMPLES THAT NEED TO BE REVIEWED IN NOTEBOOKS AND REDCAP******")
cat("\n Based on review of culture results alone, may be resolved after merge with field field.db")
cat("\n ********************************************************************************\n")
cat("\n Samples with missing sample_type or sample_type=NA:  ")
print(tbl_df(filter(culture_results, is.na(sample_type)|sample_type=="")))

cat("\n Ant Nasal samples with either passage or focus assay results (even if neg/0):  ")
print(tbl_df(culture_results%>%filter(antnasal==T, focus.i==T | passage.i==T)))

cat("\n Impactor samples with either passage or focus assay results (even if neg/0):  ")
print(tbl_df(culture_results%>%filter(impactor==T, focus.i==T | passage.i==T)))

cat("\n NP swabs with a focus assay but no passage (even invalid) assay, or with a passage but no focus assay:  ")
print(tbl_df(culture_results%>%filter(np==T, (focus.i==T & passage.i==F) | (focus.i==F & passage.i==T))))

cat("\n G-II condensate with a focus assay but no passage (even invalid) assay, or with a passage but no focus assay:  ")
print(tbl_df(culture_results%>%filter(condensate==T, (focus.i==T & passage.i==F) | (focus.i==F & passage.i==T))))

cat("\ Culture_results: \n")
print(head(tbl_df((culture_results))))
summarise(culture_results)


cat("\n\n***#7.	Merge 3***\n")

cat("\n", "***************Merge3 = merge of merge2 with (culture_results from REDCap) by sample_id (only)**************** \n")
cat("\n x=merge2, y=culture_results")
cat("\n Variables merge2:\n")
print(names(merge2))
cat("\n Variables culture_results\n")
print(names(culture_results))

merge3= merge(merge2,culture_results,c('sample_id'),all=TRUE)
cat("\n", "Number of rows in merge2: ", nrow(merge2))
cat("\n", "Number of rows in culture_results: ", nrow(culture_results))
cat("\n", "Number of rows in merge3: ", nrow(merge3))
cat("\n", "Number of cols in merge3: ", ncol(merge3), "\n")
print(head(tbl_df(merge3)))
cat("\n", "Table to check source of records after merge 3 \n")
print(ftable(addmargins(with(merge3,table(merge2.i,cr.i,exclude=c())))))

cat("\n", "Do date_visit match? \n")
d.err <- filter(merge3,date_visit.x!=date_visit.y | is.na(date_visit.x!=date_visit.y))
cat("\n Number of samples where the date_visit.x (merge2) not equal date_visit.y (culture_results)", nrow(d.err),"\n")
cat(" First 10 rows with non matching date_visit ordered by sample_id")
print(top_n(select(d.err, sample_id,date_visit.x,date_visit.y,enrolled,visit_num,sample_type.x, sample_type.y),10,sample_id))

cat("\n", "Do sample_type match? \n")
t.err <- filter(merge3, sample_type.x!=sample_type.y | is.na(sample_type.x!=sample_type.y))
cat("\n Number of rows where sample types don't match =", nrow(t.err), "\n")
cat("\n Columns show whether culture_result sample types were missing?, Rows likewise for merge2 sample type. \n")
print(ftable(addmargins(with(merge3,table(miss.x<-is.na(sample_type.x),miss.y<-is.na(sample_type.y))))))

cat("\n Table of sample types by data source. (x=merge2, y=culture_results)\n")
print(ftable(addmargins(with(merge3,table(merge2.i,cr.i,sample_type.x,sample_type.y,exclude = c())))))
cat("\n All non-matching sample types seem to be due to missing (NA) values. \n")

cat("\n", "Samples in merge3 where culture_results data had no match in merge2 (field data):\n")
x <- filter(merge3,cr.i==T,is.na(merge2.i))
print(select(x,sample_id))

cat("\n All sample ids begining with 237 in merge3 \n")
print(filter(merge3,grepl("^237",sample_id)))
cat("\n All sample ids begining with 237 in culture_results \n")
print(filter(culture_results,grepl("^237",sample_id)))
cat("\n Sample 237_6 is an enrolled roommate, enrolled based on fever, therefore second NP swab should be in the lab.\n")

cat("\n All sample ids begining with 301 in merge3\n")
print(filter(merge3,grepl("^301",sample_id)))
cat("\n All records in merge3 with field_subj_id 301 from either source dataframe\n")
print(filter(merge3,field_subj_id.x==301 | field_subj_id.y==301))
cat("\n Records for samples 301_x in culture_results: \n")
print( filter(culture_results,grepl("^301",sample_id)))



cat("\n All samples collected on same date as subject 301 \n")  
x <- select(filter(merge3, field_subj_id.y=="301"),date_visit.x)

y <- right_join(merge3,x, by="date_visit.x")
cat("\n Throat and Condensate samples from enrolled subject on the same day as 301: \n"  )
print(filter(y,sample_type.x %in% c("Throat Swab", "GII condensate NO mask") ))
cat("\n Looks like 301_3 and 301_5 are erroneous duplicative entries for 302_3 and 302_5: Will delete extra 301 samples.")
merge3 <- filter(merge3,!(sample_id %in% c("301_3","301_5")))

cat("\n", "Table to check source of records after merge 3 clean-up\n")
print(ftable(addmargins(with(merge3,table(merge2.i,cr.i,exclude=c())))))

merge3 <- mutate(merge3, 
                 subject_id = ifelse((!is.na(field_subj_id.x)|field_subj_id.x==""),field_subj_id.x,field_subj_id.y),
                 sample_type = ifelse((sample_type.x %in% c("Nasopharyngeal swab","Impactor 5 um NO mask","GII condensate NO mask",
                                                            "anterior nasal swab","Throat Swab"))
                                      ,sample_type.x,sample_type.y),
                 date_visit =  as.Date(ifelse(!is.na(date_visit.x),date_visit.x,date_visit.y),origin = "1970-01-01")
)
merge3 <- select(merge3, -contains(".x"), -contains(".y") )

cat("\n Number of rows in merge3 after clean-up with no sample type: ", nrow(filter(merge3,is.na(sample_type)|sample_type=="")),"\n")
cat("\n Number of rows in merge3 after clean-up with no subject id: ", nrow(filter(merge3,is.na(subject_id)|subject_id=="")),"\n")
cat("\n Number of rows in merge3 after clean-up with no date visit: ", nrow(filter(merge3,is.na(date_visit))),"\n")

cat("\n Samples without sample type: \n")
print(filter(merge3,is.na(sample_type)))


### - for not enrolled, figure out how to keep the one  NP swabs that was sent to the lab and discard the other record.


#   cat("\n", "Example of samples from merge2 that had no match in culture_results:\n")
#   print(head(filter(merge3,is.na(cr.i),merge2.i==T)))
#   
#   cat("\n Sample types for merge2 records that are missing in culture_results\n")
#   miss.cr <- filter(merge3,is.na(cr.i))
#   print(addmargins(with(miss.cr,(table(sample_type.x,cr.i,exclude = c())))))
#   
#   cat("\n Sample types for records that are in merge2 and culture_results\n")
#   in.both <- filter(merge3,cr.i==T,field.db1.i==T)
#   print(addmargins(with(in.both,(table(sample_type.x,sample_type.y,exclude = c())))))
#   
#   cat("\n Large summary of merge 2 \n")
#   print(ftable(addmargins(with(merge2,(table(field.db1.i,cr.i,sample_type.x,sample_type.y,exclude = c()))))))

cat("\n Examine the first visit to see which NP samples were cultured \n") 
np <- select(filter(merge3,np==T|sample_type=="Nasopharyngeal swab"), subject_id,sample_id,enrolled,visit_num,validp,ct)  
np1 <- filter(np,visit_num==1)
np1a <- select(mutate(np1,cultured=!is.na(validp)|!is.na(ct)),sample_id,subject_id,cultured)
np1a <- select(separate(np1a,sample_id,c("id","sample.id"),"_"),subject_id,sample.id,cultured)
np1a.s <- spread(np1a,sample.id,cultured)
names(np1a.s)[2:length(names(np1a.s))]<-paste("sample",names(np1a.s)[2:length(names(np1a.s))],sep="_")
print(summary(np1a.s))
cat("\n NA, means not cultured / assayed. 
    Conclude that all cultured visit 1 samples were called \"_1\" by the lab \n ")

cat("\n Examine all visits for which NP samples were cultured \n")

# -- before spread can run on all np samples, must correct for reassigning 225 to 250 to avoid duplication of sample_#
np <- mutate(np,sample_id=ifelse(subject_id==250 & visit_num==2,paste(sample_id,"a",sep=""),sample_id))
np.a <- select(mutate(np,cultured=!is.na(validp)|!is.na(ct)),sample_id,subject_id,cultured)
np.a <- select(separate(np.a,sample_id,c("id","sample.id"),"_"),subject_id,sample.id,cultured)
np.a.s <- spread(np.a,sample.id,cultured)
names(np.a.s)[2:length(names(np.a.s))]<-paste("sample",names(np.a.s)[2:length(names(np.a.s))],sep="_")
print(summary(np.a.s))
cat("Subjects with sample_6 cultured or subject_id=250 (after reassigning subject 225 to subject 250): \n")
print(filter(np.a.s,sample_6==T|subject_id==250))

cat("\n All samples for subject 247 \n")
print(filter(merge3,subject_id == "247"))
cat("\n All samples for subje_id == 250 \n")
print(filter(merge3,subject_id == "250"))

merge3 <- rename(merge3,subject.id=subject_id, sample.id=sample_id, focus.ct=ct, sample.type = sample_type, 
                 date.visit=date_visit, g2.run=g2_run, visit.num = visit_num)
samples.cc <- select(merge3,subject.id, date.visit, sample.id, sample.type, g2.run, visit.num, passpos, validp, focus.ct, enrolled,rapid_flu___3,rapid_flu_loc,body_temp)



saveRDS(samples.cc,file=paste(Out.dir,"EMIT_samples.cc.RDS",sep=""))




cat("\nProgram Merge1-3.R End runtime = ", date(),"\n")  
sink()
closeAllConnections()  
