rm(list=ls())

# R version 4.4.2 (2024-10-31)
# Platform: x86_64-pc-linux-gnu
# Running under: Ubuntu 22.04.5 LTS

library(dplyr) # dplyr_1.1.4
library(TMB) # TMB_1.9.16
library(aghq) # aghq_0.4.1
library(magrittr) # magrittr_2.0.3
library(reshape2) # reshape2_1.4.4
library(gridExtra) # gridExtra_2.3
library(cowplot) # cowplot_1.1.3
library(readr) #readr_2.1.5
library(ggplot2) #ggplot2_3.5.1
library(lubridate) #lubridate_1.9.4

compile(file="./section_3/cpp/model4.cpp")
try(dyn.unload(dynlib("./section_3/cpp/model4")),silent = TRUE)
dyn.load(dynlib("./section_3/cpp/model4"))

## Load wastewater + case counts
load("./section_4/data/cases_weekly_new_zealand.RData")
load(file="./section_4/results_model_NZ.RData")
source("./functions_general/prep_data_covid_with_fittedwastewater.R")
