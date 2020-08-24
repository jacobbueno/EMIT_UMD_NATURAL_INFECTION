# EMIT_UMD_NATURAL_INFECTION

EMIT UMD Natural Infection study. This README describes the Data cleaning and analysis project using the EMIT data from the University of Maryland campus community. The results from this study were published in PNAS (Yan et al., 2018). The goal of this R project is to take all of the raw data files from the study, clean them and merge them together into an authoritative dataset that is ready for SAS analysis with tobit mixed models. The scripts that do that are described here.

Date: December 15, 2018
Updated: February 14, 2019

The goal of this repo is to provide all of the code required to clean the raw EMIT UMD campus population data, merge the data into analyitical datasets, and complete a few brief analyses with this data. 

There are 4 cleaning/analytical .R scripts in this repo and one .R script (EMIT_UMD_Natural_Infection_Source_Scripts.R) to source and render markdown reports for these 4 cleaning/analytical .R scripts. 

There are 4 word documents which are Markdown reports that have been compiled directly from the scripts. These word documents have names that match their corresponding mother script. 

The 4 cleaning/analytical .R scripts in this repo are described below (with the script of greatest importance being the first one listed -- the EMIT_UMD_Natural_Infection_Cleaning.R script):


1) EMIT_UMD_Natural_Infection_Cleaning.R 

This is the main cleaning file for this repo. This file walks one through the entire process of selecting the various, deidentified, raw datafiles from the EMIT_UMD study, and cleaning and merging them into usable curated dataframes. 7 analytical datasets were produced from this script and they are housed in Box.com in EMIT/EMIT_Data_Analysis_Jake/EMIT_UMD_Natural_Infection/Curated Data/Analytical Datasets:

i. all_screened.csv: This is the dataset that is made by using all of the cleaned data and keeping information on all 355 of the screened study participants from the EMIT_UMD campus community (n=355 subjects; 276 g-ii sampling events)

ii. all_cases.csv: all data on all participants who met case criteria and were enrolled (n=178 subjects; 276 g-ii sampling events)

iii. all_cases_gii_samples.csv: data only where a G-II sample was taken for all participants who met case criteria and were enrolled (n=178 subjects; 276 g-ii sampling events)

iv. flu_cases.csv: all data on enrolled cases who were positive for influenza virus (n=158 subjects; 250 g-ii sampling events)

v. flu_cases_gii_samples.csv: data only where a G-II sample was taken for all enrolled cases who were positive for influenza virus (n=158 subjects; 250 g-ii sampling events)

vi. PNAS_data_full.csv: final, cleaned dataset from among enrolled cases who were positive for influenza virus and who had complete cough data, limited to visits between days 1 and 3 of illness onset, and had complete PCR data. (n=142 subjects; 218 g-ii sampling events)

vii.finaldatasetrepeatupdate.csv: (n=142 subjects; 218 g-ii sampling events) This is the dataset that is made by using the original code from Jing and Dr. Milton that creates a final dataframe that can be used in the analyses for the PNAS manuscript (these analyses were done in SAS). This dataset contains data on 142 subjects for 276 exhaled breath sampling instances. However, there appears to have been some missing script somewhere that Jing used to create her final dataframe that was used in the tobit models in SAS. The missing script adds the following variables to this final dataframe: centimeterheight, kilogramweight, BMI, cur_asthma, lung_sym_2pos, fluvac_last2y, bothyear, age. Because there was no code in the original set of scripts authored by Jing to get these variables into the final dataset, I have decided to not add them to this dataframe how the available code woud have created it, but rather, to add them in a subsequent dataframe that is another one of these 7 final analytical datasets (number vii). There are a few other minor issues detected with this original finaldatasetrepeatupdate.csv df that are described when introducing analytical dataset number vii below. 

A lucidchart diagram provides an excellent description of the cleaning and merging process that gives rise to these 7 datasets. This lucidchart diagram can be accessed here: https://www.lucidchart.com/invitations/accept/2d2c49c7-0b46-42e6-af8a-bec41ceff41f. 

Additionally, the word file under EMIT_Data_Analysis_Jake/EMIT_UMD_Natural_Infection/Analysis Notes called "samples need to be removed from final anlaysis.docx" provides a summary (see the table called "Moving from the total screened df to the PNAS analytical df" in the document) of how the dataset with data from all of the screened participants, dwindles down to the final, cleaned dataset. This is also reflected in the lucidchart diagram. 

Several key discrepancies have been raised between the resulting PNAS_data_full.csv and Jing's original PNAS cleaned dataset (which can be found in the EMIT_Data_Analysis box.com directory):

a) The labelling of the type (flu A or B) or infection for participants who were coinfected with both types, in the previous version, incorrectly labeled all of the pcr results as being either flu A or B. This was updated to reflect the accurate type of virus that was being targetting in the associated pcr result variables. This resulted in the addition of a few observations to the datasets. 

b) There were a couple of instances where pcr files had 2 different names, but were actually the same file. This occured where the date in the name of the pcr experiment used a "0" in front of a 2-digit month abbreviation (i.e., 07 for august, as opposed to simply 7). These instances were removed in the updated dataframes, thus eliminating a few, repeated observations from the datasets. 

c) The RNA copy#-to-virus particle ratio for the viral standards (PR/8 and B/Lee) used in the original cleaning process was found to be 80 and 411 for flu A and B, respectively. Michael Grantham's original experiments indicate that this ratio should actually be 250 and 272 for flu A and B, respectively. However, after applying the overall EMIT project's standard curve information to these experiments these final, adjusted ratios become 80 and 411 for flu A and B, repectively. The PCR experiments that support this are found in the "EMIT_Data_Analysis_Jake/EMIT_UMD_Natural_Infection/UMD_Raw_Data/EMIT RNA copies per virion" directory.

d) It appears that the variables cur_asthma and lung_sym_2pos in Jing's original dataset, cannot be reproduced from any of the available raw data. Email correspondence with Jing has yet to reveal how to reproduce these variables, or what they mean precisely. There is an "asthma" variable in the dataset that doesn't seem to correspond exactly with the "cur_asthma" variable in Jing's original dataset. 

e) The fluvac_last2y variable appears to have quite a few instances of missingness. NAs, as opposed to 0 (for no flu vaccine within the last 2 years), and 1 (for yes, flu vaccine taken within the last 2 years), were observed for a number of participants. Jing's original dataset had marked these NA's as 0, but I have decided to not do this unless otherwise directed. Thus, this represents another discrepancy. As a result, the bothyear variable (flu vaccine taken both years) is also incomplete because of these NA values that have not been forced to 0 as they appear to have been done in the original PNAS data.

f) Although Table S1 in Yan et al., 2018 shows that there were 178 breath collection visits from 178 subjects, we were unable to replicate this finding. We only ever show that there is breath collection data for 276 visits from 178 subjects. Perhaps there were a couple of intstances where breath collection was initiated but not completed and not marked has having occured in the REDCap database. 

****Important Note about qRT-PCR Data and RNA Copies: To be very clear on how we arrived at the final.copies variable, which shows the number of RNA copies per sample collected for NP swabs and fine and coarse aerosols, the process is as follows. The MxPro pcr assay software provides a dRn (normalized) Ct value and corresponding viral particle number, which is based on a standard curve from EM-quantified PR/8 influenza A and B/Lee virus stocks that were extracted and run on the qRT-PCR assay in the exact same way as the samples. This virus particle number is then multiplied by the RNA copy number - to - virus particle ratio that was computed from a series of experiments documented in the "EMIT_Data_Analysis_Jake/EMIT_UMD_Natural_Infection/UMD_Raw_Data/EMIT RNA copies per virion" directory. The ratios are 80 and 411 for flu A and B, respectively. Next, the RNA copy numbers were multiplied by the dilution factors based on the laboratory processing workflow (i.e., the amount of sample used for extraction and then for qRT-PCR related to the whole sample volume). For NP swabs this dilution factor was 100 (we took 50ul of sample out of 1ml total sample for extraction, and then took a fifth of the extraction volume to run the qRT-PCR assay), except for a few samples that had a dilution factor of 100 since 100ul instead of 50ul was used for the extraction. For fine and coarse aerosol samples, the dilution factor was 25 (we took 200ul of sample out of 1ml total sample for extraction, and then took a fifth of the extraction volume to run the qRT-PCR assay). Finally, these dilution-calibrated, RNA copy numbers, we multiplied by a qRT-PCR calibration factor to account for shifts in fluorescence over time (given the long period of time during which assays were run; this calibration factor was informed by standard curves and high and low interrun calibrators that were run over the course of the lab processing for these samples - details are described in the code).


2) EMIT_Enrollment_Summary.R

This script does some basic data summarizing for the study and it appears to have been useful during the ongoing participant enrollment process and for periodic checks of the data and reports. I have pulled together and cleaned up code from 2 scripts: The "Snippets analysis_1.r", and the "field database with redcap culture.R".

Of importance, the Clinical Database and the G2 Log, which are used to produce an enrollment summary which is written out as "/Users/jbueno/Box Sync/EMIT/EMIT_Data_Analysis_Jake/EMIT_UMD_Natural_Infection/Curated Data/Cleaned Data/Enrollment_Summary.csv"

This enrollment summary give s a list of the subject IDs, whether they were enrolled in the study (1 if yes, 0 if no) and a "perform_date", which is mostly likely the date of first G-II visit with each particular subject ID. 

A set of roommate pairs is also identified in this script and printed to the console, although not saved as an output dataset. This list of roommates is critical for sequence analysis to test for transmission between roommate pairs. 

Finally, a bit of data cleaning and summarizing, with comments is done with the Field Sample Database and the REDCap Database. 


3) EMIT_Enroll_Condition_for_Negative_Samples.R

This script reads in the Clinical Database (EMIT_samples.cc.RDS) and the subtyping file (EMIT_subtypes.RDS) to create a summary file that contains a list of the sampling instances with culture data (passage or quantitative focus assay) that shows that there is positive virus, despite a negative first visit NP swab (first visit NP swab is the basis of the EMIT_subtypes.RDS file). 

The output is: "/Users/jbueno/Box Sync/EMIT/EMIT_Data_Analysis_Jake/EMIT_UMD_Natural_Infection/Curated Data/Cleaned Data/negative subtype sample with positive culture or focus.RDS"

This output can be used in subsequent analyses or exploration into the dataset as it pertains to counting individuals with evidence of viral infection by both culture and qRT-PCR methods. 


4) EMIT_UMD_Natural_Infection_Analysis.R

While most of the analysis for the EMIT_UMD study was done in SAS (especially that as part of the PNAS manuscript), here are some others analyses that we will include as well. Here, we estimate the mean and variance for the ratio of ffu/copy number in fine aerosols from EMIT using a mixed model with random effect of person and fixed effect of day since onset of symptoms. 

5) metadata_prep_sequencing.R

Preparing metadata for sequencing analysis aims (CEIRS Option hypothesis 1)

