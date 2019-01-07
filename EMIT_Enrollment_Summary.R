# Title: EMIT_Enrollment_Summary
# Author: Jacob Bueno de Mesquita with material from Jing Yan and Don Milton
# Date: January 7, 2019

#### **** Using Script: Jing Yan and Dr. Milton's "Snippets analysis_1.r" **** ####

# Perhaps the earlier lines of script in this program address what "Snippets analysis_1.r" was mostly getting at.
# However, "Snippets analysis_1.r" provides some, perhaps useful, summary information and looks at roommates. 

###
# Original file information:

# By Jing Yan & Don Milton
# December 14, 2015 - December 21, 2015
# Purpose: Combine clinical redcap data with GII data from redcap and check that all subject enrolled
#          according to redcap clinical data have the appropriate number of GII records and that
#          persons not enrolled only have screening visits (up to 3) and no GII records.
#          Will also generate a list of subjects showing whether they were enrolled and how many GII
#          sessions they completed. 
#		   Will identify roommate screenings and roommates enrolled after screening as roommates. 

###

#### READ in and work with Clinical Database ####

clinical_umd <- read.csv(clinical_in_file)

# clinical_in_file was an object that was taken from reading in "EMIT_UMD_Natural_Infection/UMD_Raw_Data/REDCAP/EMITClinicalUMD2013.csv"
# This was already read in with the "Merge_1-3.R-update.r" script -- see above.

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

# g2_in_file was an object already read in earlier from 'EMIT_UMD_Natural_Infection/UMD_Raw_Data/GII/EMITGIILogUMD2013.csv'
# See above for where this was already read in with material from the "merge_1-3.R-update.R" script

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