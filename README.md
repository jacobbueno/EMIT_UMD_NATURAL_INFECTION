# EMIT_UMD_NATURAL_INFECTION

EMIT UMD Natural Infection study.

The goal of this repo is to provide all of the code required to clean the raw EMIT_UMD_campus data, merge the data into analyitical datasets, and complete a few brief analyses with this data. 

There are 4 R scripts in this repo (with the script of greatest importance being the first one listed -- the EMIT_UMD_Natural_Infection_Cleaning.R script):

1) EMIT_UMD_Natural_Infection_Cleaning.R is the main cleaning file for this repo. This file walks one through the entire process of selecting the various, deidentified, raw datafiles from the EMIT_UMD study, and cleaning and merging them into usable curated dataframes. 7 analytical datasets were produced from this script and they are housed in Box.com in EMIT/EMIT_Data_Analysis_Jake/EMIT_UMD_Natural_Infection/Curated Data/Analytical Datasets:

i. finaldatasetrepeatupdate.csv: This is the dataset that is made by using the original code from Jing and Dr. Milton that creates a final dataframe that can be used in the analyses for the PNAS manuscript (these analyses were done in SAS). This dataset contains data on 142 subjects for 276 exhaled breath sampling instances. However, there appears to have been some missing script somewhere that Jing used to create her final dataframe that was used in the tobit models in SAS. The missing script adds the following variables to this final dataframe: centimeterheight, kilogramweight, BMI, cur_asthma, lung_sym_2pos, fluvac_last2y, bothyear, age. Because there was no code in the original set of scripts authored by Jing to get these variables into the final dataset, I have decided to not add them to this dataframe how the available code woud have created it, but rather, to add them in a subsequent dataframe that is another one of these 7 final analytical datasets (number vii). There are a few other minor issues detected with this original finaldatasetrepeatupdate.csv df that are described when introducing analytical dataset number vii below. 

ii. all_screened.csv: This is the dataset that is made by using all of the cleaned data and keeping information on all 355 of the screened study participants from the EMIT_UMD campus community. 

iii. all_cases.csv: This is the dataset that is made by 

iv.

v.

vi.

vii.

2)


3)


4)




Note the adding on of the BMI and flu variables, etc. using the finaldataset_merge object like was done to the PNAS_data_full df. This can be done to the other variables. 

Also note the issue with the vaccination data and how we believe that some of it was missing. 






