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

sessionInfo()


