rm(list=ls())

# R version 4.4.2 (2024-10-31)
# Platform: x86_64-pc-linux-gnu
# Running under: Ubuntu 22.04.5 LTS

library(dplyr) # dplyr_1.1.4
library(TMB) # TMB_1.9.16
library(magrittr) # magrittr_2.0.3
library(reshape2) # reshape2_1.4.4
library(gridExtra) # gridExtra_2.3
library(cowplot) # cowplot_1.1.3
library(lubridate) # lubridate_1.9.4
library(ggplot2) # ggplot2_3.5.1
library(zoo) #zoo_1.8-12
library(readr) # readr_2.1.5

compile(file="./section_3/cpp/model4.cpp")
try(dyn.unload(dynlib("./section_3/cpp/model4")),silent = TRUE)
dyn.load(dynlib("./section_3/cpp/model4"))

## Load wastewater + case counts
load("./section_4/data/cases_weekly_new_zealand.RData")
load(file="./section_4/results_model_NZ.RData")
source("./functions_general/prep_data_covid_with_fittedwastewater.R")
source('./functions_general/process_results_epidemic.R')


sites <- read_csv("./section_4/data/sites.csv")
sites <- sites %>% 
  dplyr::select(SampleLocation, Region, Population) %>% 
  rename("site_id"= SampleLocation,
         "weight"=Population)

sites2 <- sites %>% 
  group_by(Region) %>% 
  mutate(region_weight = sum(weight)) %>% 
  mutate(region_weight = ifelse(Region == "Auckland", region_weight-weight[which(site_id == "AU_Mangere")], region_weight)) %>% 
  mutate(weight = ifelse(site_id == "AU_Mangere", 1.716*1000000- region_weight, weight))

data_foranalysis <- prep_data(case_data = cases_epiweekly,
                              y_var = "number_of_cases",
                              results = results,
                              AR = TRUE,
                              weight_ratio = FALSE,
                              weights_datadic = sites2,
                              date_forfilter_strict_upr = ymd("2023-03-31"))


### Fit model 
MI_models <- list(NULL)

set.seed(2917)
ss <- sample(1:3000,500)
for (i in 1:length(ss)){
  
  tmbdat <- data_foranalysis$tmbdat[[ss[i]]]
  tmbdat$case_counts <- c(tmbdat$obs_start_case, tmbdat$case_counts)
  tmbdat$ratio <- c(1, tmbdat$ratio_v_u_fixed)
  tmbdat$lag = 0
  tmbdat$u1 = 1
  tmbdat$alpha1 = 0.5 #median
  tmbdat$mean_z0 = tmbdat$obs_start_case
  tmbdat$sd_z0 = 500
  
  tmbparams <- list(
    W = c(rep(0, length(tmbdat$case_counts)-tmbdat$lag), 
          rep(400,length(tmbdat$case_counts)-tmbdat$lag)), # W = c(U,beta,Z); 
    theta_p=1
  )
  
  
  ff <- TMB::MakeADFun(
    data = tmbdat,
    parameters = tmbparams,
    random = "W",
    DLL = "model4",
    silent = TRUE
  )
  
  aghq_k = 10
  
  mdl1 <- aghq::marginal_laplace_tmb(ff,k=aghq_k,startingvalue = c(1))
  samps1 <- aghq::sample_marginal(mdl1, M = 3000) 
  
  MI_models[[i]] <- samps1
  print(i)
}

save(file="./section_4/results_popweighted/results_NZ_v_u_fixed.RData", list = "MI_models")


load("./section_4/results_popweighted/results_NZ_v_u_fixed.RData")
results <- process_results(data_foranalysis = data_foranalysis,
                           MI_models = MI_models,
                           full_timeseries = FALSE,
                           adj = FALSE)

save(file="./section_4/results_popweighted/results_processed_NZ_v_u_fixed.RData", list = "results")


