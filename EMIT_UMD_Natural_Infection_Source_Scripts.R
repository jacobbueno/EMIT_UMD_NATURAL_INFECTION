# Title: EMIT_UMD_Natural_Infection_Source_Scripts.R
# Author: Jacob Bueno de Mesquita
# Date: February 7, 2019
# Objective: Run each of the 4 scripts required to reproduce the EMIT UMD Natural Infection analyses and create markdown reports for each of them. Remember that the main objective of the scripts in this project is to explore the data and produce clean datasets that can be imported into SAS for tobit model with mixed effects. Tobit model with mixed effects cannot be done in R (until someone makes a new package or updates one that has that capacity).

library(markdown)
library(rmarkdown)
library(htmlTable)

source("/Users/jbueno/Desktop/Dissertation_git/EMIT_UMD_NATURAL_INFECTION/EMIT_UMD_Natural_Infection_Cleaning.R")
render("/Users/jbueno/Desktop/Dissertation_git/EMIT_UMD_NATURAL_INFECTION/EMIT_UMD_Natural_Infection_Cleaning.R", "word_document")

source("/Users/jbueno/Desktop/Dissertation_git/EMIT_UMD_NATURAL_INFECTION/EMIT_Enrollment_Summary.R")
render("/Users/jbueno/Desktop/Dissertation_git/EMIT_UMD_NATURAL_INFECTION/EMIT_Enrollment_Summary.R", "word_document")

source("/Users/jbueno/Desktop/Dissertation_git/EMIT_UMD_NATURAL_INFECTION/EMIT_Enroll_Condition_for_Negative_Samples.R")
render("/Users/jbueno/Desktop/Dissertation_git/EMIT_UMD_NATURAL_INFECTION/EMIT_Enroll_Condition_for_Negative_Samples.R", "word_document")

source("/Users/jbueno/Desktop/Dissertation_git/EMIT_UMD_NATURAL_INFECTION/UMD_Natural_Infection_Analysis.R")
render("/Users/jbueno/Desktop/Dissertation_git/EMIT_UMD_NATURAL_INFECTION/UMD_Natural_Infection_Analysis.R", "word_document")
