rm(list=ls())

# R version 4.4.2 (2024-10-31)
# Platform: x86_64-pc-linux-gnu
# Running under: Ubuntu 22.04.5 LTS

library(dplyr) # dplyr_1.1.4
library(TMB) # TMB_1.9.16
library(aghq) # aghq_0.4.1
library(magrittr) # magrittr_2.0.3
library(reshape2) # reshape2_1.4.4

load("./section_3/data/work_d_toronto_2024_08_08.RData")
load(file="./section_3/results_model_ON_phachist.RData")
source("./functions_general/prep_data_covid_with_fittedwastewater.R")

compile(file="./section_3/cpp/model4.cpp")
try(dyn.unload(dynlib("./section_3/cpp/model4")),silent = TRUE)
dyn.load(dynlib("./section_3/cpp/model4"))


weights_datadic = data.frame(site_id = c("TAB","THC","THU","TNT"), 
                             weight = c(0.499/(0.177+0.232+0.064+0.499),
                                        0.177/(0.177+0.232+0.064+0.499),
                                        0.232/(0.177+0.232+0.064+0.499),
                                        0.064/(0.177+0.232+0.064+0.499)))

data_foranalysis_full <- prep_data(case_data = work_d_toronto,
                                   y_var = "number_of_cases",
                                   results = results,
                                   AR=TRUE,
                                   weight_ratio = TRUE,
                                   weights_datadic = weights_datadic)



m1 = data_foranalysis_full$analysis_d[[1]] %>% 
  filter(earliest_week_end_date <= "2021-12-01") %$% earliest_week_end_date %>% max()

m2 = data_foranalysis_full$analysis_d[[1]]$earliest_week_end_date %>% max()

seq(m1, m2, 7) %>% length()

n1 = data_foranalysis_full$analysis_d[[1]] %>% 
  filter(earliest_week_end_date <= m1) %$% earliest_week_end_date %>% length()


#######
for (j in seq(0,126,6)){
  data_foranalysis <- as.list(NULL)  
  data_foranalysis$tmbdat <-  lapply(data_foranalysis_full$tmbdat, 
                                     function(x){
                                       x$case_counts <- x$case_counts[1:(n1 + j)]
                                       x$ratio <- x$ratio[1:(n1 + j)]
                                       x$ratio_v_u <- x$ratio_v_u[1:(n1 + j)]
                                       x$ratio_v_u_fixed <- x$ratio_v_u_fixed[1:(n1 + j)]
                                       return(x)
                                     })
  
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
            rep(2000,length(tmbdat$case_counts)-tmbdat$lag)), # W = c(U,beta,Z); 
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
    samps1 <- aghq::sample_marginal(mdl1, M = 500) 
    
    MI_models[[i]] <- samps1
    print(i)
  }
  
  save(file=paste0("./section_3/results_pred_popweighted_fixed/results_",j,".RData"), list = "MI_models")
}

#######

#######
results_pred <- NULL
for (j in seq(0,126,6)){
  pred_week = 5
  load(paste0("./section_3/results_pred_popweighted_fixed/results_",j,".RData"))
  
  ncounts = n1 + j
  lag = 0
  
  set.seed(29172)
  len = ncounts+1 - lag 
  logit_p_samps1 <- lapply(MI_models, function(samps){samps$samps[1:(len),]})
  Z_samps1 <- lapply(MI_models, function(samps){samps$samps[(len  + 1):((len  +len)),]})
  thetasamps <- lapply(MI_models, function(samps){samps$thetasamples[[1]]})
  
  lensim = ncol(Z_samps1[[1]])*length(Z_samps1)
  
  logit_p_samps1 <- Reduce(cbind, logit_p_samps1)
  thetasamps <- Reduce(c, thetasamps)
  
  set.seed(2)
  logit_p_sampspred <- matrix(data = logit_p_samps1[nrow(logit_p_samps1),], nrow = 1)
  for (i in 2:(pred_week+1)){
    logit_p_sampspred <- rbind(logit_p_sampspred,
                               rnorm(ncol(logit_p_sampspred), logit_p_sampspred[(i-1),], sd = exp(-0.5*thetasamps)))}
  
  p_sampspred = exp(logit_p_sampspred)/(1+  exp(logit_p_sampspred))
  
  
  Z_samps1 <- Reduce(cbind, Z_samps1)
  set.seed(2917)
  ss <- sample(1:3000,500)
  alphas = t(Reduce(rbind,lapply(data_foranalysis_full$tmbdat[ss], function(x) matrix(rep(x$ratio_v_u_fixed[(n1 + j + 1):(n1+j + pred_week)],500), byrow=TRUE, ncol = pred_week ))))
  
  set.seed(2)
  Z_sampspred <- matrix(data = Z_samps1[nrow(Z_samps1),], nrow = 1)
  Z_sampspred <- abs(Z_sampspred)
  
  for (i in 2:(pred_week+1)){
    Z_sampspred <- rbind(Z_sampspred,
                         rnorm(ncol(Z_sampspred), alphas[i-1,]*Z_sampspred[(i-1),], sd = sqrt(alphas[i-1,]*Z_sampspred[(i-1),]) ))
    Z_sampspred <- abs(Z_sampspred)}
  
  set.seed(3)
  y_pred = matrix(NA, nrow = pred_week, ncol = lensim)
  for (i in 1:pred_week){
    y_pred[i, ] = rnorm(lensim, p_sampspred[(i+1),]*Z_sampspred[(i+1),], sqrt(p_sampspred[(i+1),]*(1-p_sampspred[(i+1),])*Z_sampspred[(i+1),]))
  }
  
  pred_df =  data_foranalysis_full$analysis_d[[1]] %>% 
    dplyr::select(week_start_date, earliest_week_end_date, y, total_number_of_tests, admissions, deaths) %>% 
    slice((n1 + j + 1):(n1+j + pred_week)) 
  
  pred_df$pred_y_med = as.numeric(apply(y_pred, MARGIN=1,median))
  pred_df$pred_y_95_upr = as.numeric(apply(y_pred, MARGIN=1,quantile,0.975))
  pred_df$pred_y_95_lwr = as.numeric(apply(y_pred, MARGIN=1,quantile,0.025))
  pred_df$pred_y_80_upr = as.numeric(apply(y_pred, MARGIN=1,quantile,0.9))
  pred_df$pred_y_80_lwr = as.numeric(apply(y_pred, MARGIN=1,quantile,0.1))
  
  pred_df$pred_prob_med = as.numeric(apply(p_sampspred[-1,], MARGIN=1,median))
  pred_df$pred_prob_95_upr = as.numeric(apply(p_sampspred[-1,], MARGIN=1,quantile,0.975))
  pred_df$pred_prob_95_lwr = as.numeric(apply(p_sampspred[-1,], MARGIN=1,quantile,0.025))
  pred_df$pred_prob_80_upr = as.numeric(apply(p_sampspred[-1,], MARGIN=1,quantile,0.9))
  pred_df$pred_prob_80_lwr = as.numeric(apply(p_sampspred[-1,], MARGIN=1,quantile,0.1))
  
  
  pred_df$pred_Z_med = as.numeric(apply(Z_sampspred[-1,], MARGIN=1,median))
  pred_df$pred_Z_95_upr = as.numeric(apply(Z_sampspred[-1,], MARGIN=1,quantile,0.975))
  pred_df$pred_Z_95_lwr = as.numeric(apply(Z_sampspred[-1,], MARGIN=1,quantile,0.025))
  pred_df$pred_Z_80_upr = as.numeric(apply(Z_sampspred[-1,], MARGIN=1,quantile,0.9))
  pred_df$pred_Z_80_lwr = as.numeric(apply(Z_sampspred[-1,], MARGIN=1,quantile,0.1))
  
  pred_df$model = paste0("model",j)
  
  results_pred <- rbind(results_pred, pred_df)
}


results_pred %>%
  ggplot(aes(week_start_date, y))+
  geom_point()+
  geom_point(aes(week_start_date, pred_y_med, col = model), show.legend = FALSE)+
  geom_errorbar(aes(week_start_date, ymin = pred_y_95_lwr, ymax = pred_y_95_upr, col = model), show.legend = FALSE)+
  scale_y_continuous(trans = "log10")

save(file="./section_3/results_pred_popweighted_fixed/results_pred.RData", results_pred)
